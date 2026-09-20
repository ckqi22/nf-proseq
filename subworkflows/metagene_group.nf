#!/usr/bin/env nextflow
//
// SUBWORKFLOW: metagene_group
// deepTools 版 TSS metagene（按组出图版，与按样本版 metagene.nf 并列）：
//   - BIGWIGAVERAGE: 组内多样本 bigwigAverage 逐碱基平均（重复合并）
//   - PLOTPROFILE: 全基因平均 TSS metagene 曲线（仅 sense；每组一张图、一条线）
//   - PLOTHEATMAP: 逐基因 TSS 热图（每行一个基因，按信号升序；每组一张）
//
// 输入通道 bigwig_cpm 与按样本版相同（每样本一个 tuple：meta + 5'-端 CPM bigWig 对），
// meta.group 指定分组。仅「≥2 样本的显式分组」出组图：'unknown'（samplesheet 缺
// group 列的缺省值）与单样本组在 groupTuple 后被过滤，不产生组图。
// 组内样本先平均成每组一对 bigWig，再走 COMPUTEMATRIX（正/负链基因分开
// computeMatrix，minus 传 --scale -1 翻正负值）→ computeMatrixOperations rbind 单 group 矩阵。
//

// DSL2 同一 process 不能在单个 workflow 中调用两次 → 用别名各调一次
include { BIGWIGAVERAGE                           } from '../modules/deeptools/bigwigAverage'
include { COMPUTEMATRIX as COMPUTEMATRIX_PROFILE  } from '../modules/deeptools/computeMatrix'
include { COMPUTEMATRIX as COMPUTEMATRIX_HEATMAP  } from '../modules/deeptools/computeMatrix'
include { COMPUTEMATRIX as COMPUTEMATRIX_COMBINED } from '../modules/deeptools/computeMatrix'
include { PLOTPROFILE                             } from '../modules/deeptools/plotProfile'
include { PLOTPROFILE as PLOTPROFILE_COMBINED     } from '../modules/deeptools/plotProfile'
include { PLOTHEATMAP                             } from '../modules/deeptools/plotHeatmap'

workflow metagene_group {
    take:
    bigwig_cpm   // channel: tuple(meta, plus_bw, minus_bw) × 样本数 — 5'-端 CPM bigWigs
    tss_bed      // channel: path(tss.bed)（BED6，第 6 列 strand）
    chrom_sizes  // channel: path(chrom.sizes)（主基因组，供 BIGWIGAVERAGE 写回 bigWig）

    main:
    def pu = params.metagene.profile_upstream   ?: 1000   // profile TSS 上游（bp）
    def pd = params.metagene.profile_downstream ?: 1000   // profile TSS 下游（bp）
    def hu = params.metagene.heatmap_upstream   ?: 50     // heatmap 上游
    def hd = params.metagene.heatmap_downstream ?: 150    // heatmap 下游
    def bs = params.metagene.bin_size           ?: 10     // bin 大小（bp）

    // 组键：'unknown'（samplesheet 缺 group 列的缺省值，main.nf:71 写入）= 未分组
    // → 退化为按样本各自成组，随后被「≥2 样本」过滤 → 不产生组图
    // 从输入 bigwig 文件名派生信号描述符 sig（无链向）与输出名后缀：
    //   S1_single_plus_cpm → plus_suffix=full 不对，这里 single_plus_cpm；sig=single_cpm；
    //   S1_full_plus_spike → plus_suffix=full_plus_spike；sig=full_spike。
    // plus_suffix/minus_suffix 供 BIGWIGAVERAGE 拼组平均输出名，sig 供下游
    //   computeMatrix/plot 输出命名 + y 轴标签（与按样本版 metagene.nf 一致）。
    ch_bw = bigwig_cpm.map { meta, p, m ->
        def g = (meta.group && meta.group != 'unknown') ? meta.group : meta.sample
        def plus_suffix  = p.baseName.substring(meta.sample.size() + 1)
        def minus_suffix = m.baseName.substring(meta.sample.size() + 1)
        def sig = plus_suffix.replace('_plus', '')
        [[sample: meta.sample, group: g, sig: sig, plus_suffix: plus_suffix, minus_suffix: minus_suffix], p, m]
    }

    // 按组收集 → 过滤单样本组 → bigwigAverage 组内平均（重复合并）→ 每组一对平均 bigWig
    // 组内同信号类型 → 取 metas[0] 的 suffix/sig 即可（一致）。
    ch_avg = ch_bw.map { meta, p, m -> [meta.group, meta, p, m] }
                  .groupTuple(by: 0)
                  .filter { g, metas, plusses, minusses -> metas.size() > 1 }
                  .map { g, metas, plusses, minusses ->
                      [[sample: g, group: g, sig: metas[0].sig,
                        plus_suffix: metas[0].plus_suffix, minus_suffix: metas[0].minus_suffix],
                       plusses, minusses]
                  }
    avg_bw = BIGWIGAVERAGE(ch_avg, chrom_sizes)

    // ---- profile 矩阵（上游 pu / 下游 pd）----
    profile_matrix = COMPUTEMATRIX_PROFILE(avg_bw, tss_bed, pu, pd, bs, 'profile')

    // ---- heatmap 矩阵（上游 hu / 下游 hd）----
    heatmap_matrix = COMPUTEMATRIX_HEATMAP(avg_bw, tss_bed, hu, hd, bs, 'heatmap')

    // COMPUTEMATRIX 只有一个 emit → 调用结果即输出通道本身（单输出进程无 .out）
    PLOTPROFILE(profile_matrix)
    PLOTHEATMAP(heatmap_matrix)

    // ---- 叠加图：所有组一张 profile（组平均 bigWig 多样本矩阵 → 一图多线）----
    // 收集全部组平均 bigWig → 单个 List 传给 COMPUTEMATRIX；无 ≥2 样本分组时 avg_bw
    // 为空 → toList 得空列表，filter 跳过（不产空 all_groups 图）
    combine_all_groups = avg_bw
        .map { meta, p, m -> [meta.sample, meta.sig, p, m] }
        .toList()
        .filter { l -> l.size() > 0 }
        .map { list ->
            def sorted = list.sort { it[0] }   // 按组/样本名稳定排序 → 列表顺序确定 → -resume 可缓存
            [[sample: 'all_groups', plot_type: 'lines', per_group: true, sig: sorted[0][1], sample_names: sorted.collect { x -> x[0] }],
             sorted.collect { x -> x[2] }, sorted.collect { x -> x[3] }]
        }
    cm_combined = COMPUTEMATRIX_COMBINED(combine_all_groups, tss_bed, pu, pd, bs, 'profile')
    PLOTPROFILE_COMBINED(cm_combined)

    emit:
    // 裸 path（无 meta 包装），与按样本版一致 → main.nf 输出块（09.TSS_Metagene/）
    // 直接消费。merged 原始矩阵与 heatmap sorted_regions.bed 并入 matrix 一并交付。
    avg_bigwig = avg_bw   // tuple(meta, ${group}${plus_suffix}.bigWig, ${group}${minus_suffix}.bigWig) — 组平均 bigWig（后缀随信号）
    plot   = PLOTPROFILE.out.profile.mix(PLOTHEATMAP.out.plot)
                  .map { _meta, f -> f }
    matrix = PLOTPROFILE.out.matrix
                  .mix(PLOTHEATMAP.out.matrix,
                       PLOTHEATMAP.out.sorted_regions,
                       profile_matrix,
                       heatmap_matrix)
                  .map { _meta, f -> f }
    combined_plot   = PLOTPROFILE_COMBINED.out.profile.map { _meta, f -> f }
    combined_matrix = PLOTPROFILE_COMBINED.out.matrix.map  { _meta, f -> f }
}
