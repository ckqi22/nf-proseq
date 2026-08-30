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

    // 同一份逐碱基 bedGraph 供两个区域计数进程消费，并 emit 出去（Nextflow 多消费者广播）。
    // 区域名（type）已随 region BED 以 tuple 形式传入，不再单独硬编码。
    SINGLEBASE_COUNT_PROMOTER(GENOMECOV.out.bedgraph, promoter_bed)
    SINGLEBASE_COUNT_GENEBODY(GENOMECOV.out.bedgraph, genebody_bed)

    MERGE_PROMOTER(SINGLEBASE_COUNT_PROMOTER.out.counts.collect(), 'promoter', 'pol2_promoter')
    MERGE_GENEBODY(SINGLEBASE_COUNT_GENEBODY.out.counts.collect(), 'genebody', 'pol2_genebody')

    emit:
    bedGraph  = GENOMECOV.out.bedgraph                      // tuple(meta, plus_bg, minus_bg)  逐碱基信号
    bigwig    = GENOMECOV.out.bigwig                        // tuple(meta, plus_bw, minus_bw)  gene-strand
    promoter_counts = SINGLEBASE_COUNT_PROMOTER.out.counts  // per-sample single-base promoter counts
    promoter_matrix = MERGE_PROMOTER.out.matrix             // single-base promoter matrix
    genebody_counts = SINGLEBASE_COUNT_GENEBODY.out.counts  // per-sample single-base genebody counts
    genebody_matrix = MERGE_GENEBODY.out.matrix             // single-base genebody matrix (PI)
}
