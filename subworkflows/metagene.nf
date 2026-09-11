#!/usr/bin/env nextflow
//
// SUBWORKFLOW: metagene
// deepTools 版 TSS metagene（替代 R 版 metagene.nf）：
//   - PLOTPROFILE: 全基因平均 TSS metagene 曲线（仅 sense，单条线）
//   - PLOTHEATMAP: 逐基因 TSS 热图（每行一个基因，按信号升序）
//
// 与 R 版输入完全相同（bigwig_cpm 5'-端 CPM bigWig 对 + tss.bed）。
// COMPUTEMATRIX 内部：正/负链基因分开 computeMatrix（minus 传 --scale -1
// 翻正负值），再 computeMatrixOperations rbind 合并为单 group 矩阵。
//

// DSL2 同一 process 不能在单个 workflow 中调用两次 → 用别名各调一次
include { COMPUTEMATRIX as COMPUTEMATRIX_PROFILE } from '../modules/deeptools/computeMatrix'
include { COMPUTEMATRIX as COMPUTEMATRIX_HEATMAP } from '../modules/deeptools/computeMatrix'
include { PLOTPROFILE }   from '../modules/deeptools/plotProfile'
include { PLOTHEATMAP }   from '../modules/deeptools/plotHeatmap'

workflow metagene {
    take:
    bigwig_cpm   // channel: tuple(meta, plus_bw, minus_bw) — 5'-端 CPM bigWigs
    tss_bed      // channel: path(tss.bed)（BED6，第 6 列 strand）

    main:
    def pw = params.metagene.tss_window         ?: 1000   // profile 窗口（±pw）
    def hu = params.metagene.heatmap_upstream   ?: 50     // heatmap 上游
    def hd = params.metagene.heatmap_downstream ?: 150    // heatmap 下游
    def bs = params.metagene.bin_size           ?: 10     // bin 大小（bp）

    // ---- profile 矩阵（±pw）----
    profile_matrix = COMPUTEMATRIX_PROFILE(bigwig_cpm, tss_bed, pw, pw, bs, 'profile')

    // ---- heatmap 矩阵（上游 hu / 下游 hd）----
    heatmap_matrix = COMPUTEMATRIX_HEATMAP(bigwig_cpm, tss_bed, hu, hd, bs, 'heatmap')

    // COMPUTEMATRIX 只有一个 emit → 调用结果即输出通道本身（单输出进程无 .out）
    PLOTPROFILE(profile_matrix)
    PLOTHEATMAP(heatmap_matrix)

    emit:
    // 裸 path（无 meta 包装），与 R 版 metagene.nf 的 emit 契约一致 → main.nf
    // 输出块（09.TSS_Metagene/）直接消费。merged 原始矩阵（2 条
    // *_TSS_merged_matrix.gz）与 heatmap sorted_regions.bed 并入 matrix 一并交付。
    plot   = PLOTPROFILE.out.profile.mix(PLOTHEATMAP.out.plot)
                  .map { _meta, f -> f }
    matrix = PLOTPROFILE.out.matrix
                  .mix(PLOTHEATMAP.out.matrix,
                       PLOTHEATMAP.out.sorted_regions,
                       profile_matrix,
                       heatmap_matrix)
                  .map { _meta, f -> f }
}
