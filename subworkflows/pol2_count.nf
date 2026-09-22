#!/usr/bin/env nextflow
//
// SUBWORKFLOW: pol2_count
// Pol II 活性位点覆盖(bedtools genomecov，gene-strand 约定)：
//   GENOMECOV 从 r1_bam(read1 单端化 BAM)生成 bedGraph + bigWig，按 signal_mode 门控：
//     mode=single → 只跑单碱基活性位点端(_single；末端由 params.strandedness 决定：reverse→5'端 / forward→3'端)；
//     mode=full → 只跑全长覆盖度(_full)；
//     mode=both → 两个 alias 都跑(开发者对比 full-read vs 单碱基)。
//   复用同一份 GENOMECOV bedGraph：
//     - 逐碱基信号(bedGraph + bigWig) → 交付 / tss_meta
//     - promoter.bed / genebody.bed 单碱基区域计数 → PI
//

include { GENOMECOV                                   } from '../modules/bedtools/genomecov.nf'
include { GENOMECOV as GENOMECOV_FULL                 } from '../modules/bedtools/genomecov.nf'
include { SINGLEBASE_COUNT as SINGLEBASE_COUNT_PROMOTER } from '../modules/singlebase_count.nf'
include { SINGLEBASE_COUNT as SINGLEBASE_COUNT_GENEBODY } from '../modules/singlebase_count.nf'
include { SINGLEBASE_COUNT_MERGE as MERGE_PROMOTER      } from '../modules/singlebase_count_merge.nf'
include { SINGLEBASE_COUNT_MERGE as MERGE_GENEBODY      } from '../modules/singlebase_count_merge.nf'

workflow pol2_count {
    take:
    r1_bam        // channel: tuple(meta, r1.bam) — read1 单端化 BAM（align_bowtie2 产出）
    promoter_bed  // channel: tuple val(name), path(promoter.bed) — 最长 transcript TSS 窗口
    genebody_bed  // channel: tuple val(name), path(genebody.bed) — 最长 transcript gene body

    main:
    // signal_type 门控：mode 保留原始 signal_mode（single|full|both），use_full 决定 PI/下游取哪套。
    def mode = params.signal_mode?.trim() ?: 'single'
    def use_full = (mode == 'full')

    // 归一化 scale（对齐 normalize.R）：cpm 恒 = 1e6/total_mapped；spike = 1e6/spike_count（仅 spike 开时）。
    //   cpm/spike 均无 length 项（length 留给未来 fpkm/rpkm 扩展，当前不触发）。
    r1_bam_scaled = r1_bam.map { meta, bam ->
        def mm = meta.main_mapped ?: 0
        if (mm <= 0) { error "main_mapped <= 0 for ${meta.sample}" }
        def scale_cpm   = String.format('%.12g', 1e6d / mm)
        def scale_spike = (meta.spike_count != null && meta.spike_count > 0)
            ? String.format('%.12g', 1e6d / meta.spike_count)
            : null
        [meta + [scale_cpm: scale_cpm, scale_spike: scale_spike], bam]
    }

    // 按 mode 门控：不需要的信号喂空通道 → 进程仍被调用（.out 恒有定义、不报 "not been invoked"），
    // 空输入不产生 task（零成本）；emit / PI 直接引用 .out，未调信号即空通道。
    GENOMECOV(     (mode != 'full')   ? r1_bam_scaled : channel.empty(), 'single')
    GENOMECOV_FULL((mode != 'single') ? r1_bam_scaled : channel.empty(), 'full')

    // PI 计数用 raw bigWig（singlebase_count 里 bigWigAverageOverBed 求和 + int 取整，
    // 需整数计数，不能用 CPM/spike 浮点版）
    SINGLEBASE_COUNT_PROMOTER(use_full ? GENOMECOV_FULL.out.bigwig : GENOMECOV.out.bigwig, promoter_bed)
    SINGLEBASE_COUNT_GENEBODY(use_full ? GENOMECOV_FULL.out.bigwig : GENOMECOV.out.bigwig, genebody_bed)

    promoter_counts_ch = SINGLEBASE_COUNT_PROMOTER.out.counts
    genebody_counts_ch = SINGLEBASE_COUNT_GENEBODY.out.counts

    MERGE_PROMOTER(promoter_counts_ch.map { _meta, f -> f }.collect(),
                   promoter_counts_ch.map { meta, _f -> meta.sample }.collect(),
                   'pol2_promoter')
    MERGE_GENEBODY(genebody_counts_ch.map { _meta, f -> f }.collect(),
                   genebody_counts_ch.map { meta, _f -> meta.sample }.collect(),
                   'pol2_genebody')

    emit:
    bedGraph            = GENOMECOV.out.bedgraph                // tuple(meta, _single_plus.bedgraph, _single_minus.bedgraph)
    bigwig              = GENOMECOV.out.bigwig                  // 单碱基 raw
    bedgraph_cpm        = GENOMECOV.out.bedgraph_cpm            // 单碱基 cpm bedgraph
    bigwig_cpm          = GENOMECOV.out.bigwig_cpm              // 单碱基 cpm
    bigwig_spike        = GENOMECOV.out.bigwig_spike            // 单碱基 spike（开 spike 才有）
    bigwig_full         = GENOMECOV_FULL.out.bigwig             // 全长 raw
    bigwig_full_cpm     = GENOMECOV_FULL.out.bigwig_cpm         // 全长 cpm（取代旧 bigwig_coverage_cpm）
    bigwig_full_spike   = GENOMECOV_FULL.out.bigwig_spike       // 全长 spike（开 spike 才有）
    promoter_counts     = SINGLEBASE_COUNT_PROMOTER.out.counts  // per-sample single-base promoter counts
    promoter_matrix     = MERGE_PROMOTER.out.matrix             // single-base promoter matrix
    genebody_counts     = SINGLEBASE_COUNT_GENEBODY.out.counts  // per-sample single-base genebody counts
    genebody_matrix     = MERGE_GENEBODY.out.matrix             // single-base genebody matrix
}
