#!/usr/bin/env nextflow
//
// SUBWORKFLOW: pol2_count
// Pol II 单碱基 5' 端活性位点覆盖(bedtools genomecov，gene-strand 约定):
//   GENOMECOV 从 BAM 生成单碱基 bedGraph + bigWig：
//     _plus = + 链基因信号（正）；_minus = - 链基因信号（取负）
//
// 复用同一份 GENOMECOV bedGraph：
//   - 逐碱基信号（bedGraph + bigWig）→ 交付 / tss_meta
//   - promoter.bed / genebody.bed 单碱基区域计数 → PI
//

include { GENOMECOV                                     } from '../modules/bedtools/genomecov.nf'
include { SINGLEBASE_COUNT as SINGLEBASE_COUNT_PROMOTER } from '../modules/singlebase_count.nf'
include { SINGLEBASE_COUNT as SINGLEBASE_COUNT_GENEBODY } from '../modules/singlebase_count.nf'
include { SINGLEBASE_COUNT_MERGE as MERGE_PROMOTER      } from '../modules/singlebase_count_merge.nf'
include { SINGLEBASE_COUNT_MERGE as MERGE_GENEBODY      } from '../modules/singlebase_count_merge.nf'

workflow pol2_count {
    take:
    bam           // channel: tuple(meta, bam)
    promoter_bed  // channel: tuple val(name), path(promoter.bed) — 最长 transcript TSS 窗口
    genebody_bed  // channel: tuple val(name), path(genebody.bed) — 最长 transcript gene body

    main:
    GENOMECOV(bam)

    // signal_mode：single=单碱基5'端(默认) / full=full read 全长覆盖度
    // PI 计数用 raw bigWig（singlebase_count 里 bigWigAverageOverBed 求和 + int 取整，
    // 需整数计数，不能用 CPM 浮点版）
    def mode = (params.signal_mode?.trim() ?: 'single') == 'full' ? 'full' : 'single'

    SINGLEBASE_COUNT_PROMOTER(mode == 'full' ? GENOMECOV.out.bigwig_full : GENOMECOV.out.bigwig, promoter_bed)
    SINGLEBASE_COUNT_GENEBODY(mode == 'full' ? GENOMECOV.out.bigwig_full : GENOMECOV.out.bigwig, genebody_bed)

    promoter_counts_ch = SINGLEBASE_COUNT_PROMOTER.out.counts
    genebody_counts_ch = SINGLEBASE_COUNT_GENEBODY.out.counts

    MERGE_PROMOTER(promoter_counts_ch.map { _meta, f -> f }.collect(),
                   promoter_counts_ch.map { meta, _f -> meta.sample }.collect(),
                   'pol2_promoter')
    MERGE_GENEBODY(genebody_counts_ch.map { _meta, f -> f }.collect(),
                   genebody_counts_ch.map { meta, _f -> meta.sample }.collect(),
                   'pol2_genebody')

    emit:
    bedGraph  = GENOMECOV.out.bedgraph                      // tuple(meta, plus_bg, minus_bg)  逐碱基信号（原始）
    bigwig    = GENOMECOV.out.bigwig                        // tuple(meta, plus_bw, minus_bw)  gene-strand（原始）
    bedgraph_cpm = GENOMECOV.out.bedgraph_cpm               // tuple(meta, plus_cpm.bedgraph, minus_cpm.bedgraph)  标准化
    bigwig_cpm   = GENOMECOV.out.bigwig_cpm                 // tuple(meta, plus_cpm.bigWig, minus_cpm.bigWig)      标准化
    bigwig_coverage_cpm = GENOMECOV.out.bigwig_coverage_cpm // tuple(meta, forward_cpm.bigWig, reverse_cpm.bigWig)  覆盖度（配文献 GSE264740）
    promoter_counts = SINGLEBASE_COUNT_PROMOTER.out.counts  // per-sample single-base promoter counts
    promoter_matrix = MERGE_PROMOTER.out.matrix             // single-base promoter matrix
    genebody_counts = SINGLEBASE_COUNT_GENEBODY.out.counts  // per-sample single-base genebody counts
    genebody_matrix = MERGE_GENEBODY.out.matrix             // single-base genebody matrix
}
