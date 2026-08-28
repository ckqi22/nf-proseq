#!/usr/bin/env nextflow
//
// SUBWORKFLOW: quantification
// Per-sample featureCounts over two regions per gene (promoter pause / gene body),
// then merge each region's per-sample count files into one clean count matrix.
//
// featureCounts runs once per sample per region — the -p flag is decided per
// sample, so single-end and paired-end samples coexist correctly in one run.
//

include { EXTRACT_R1                              } from '../modules/samtools/extract_r1.nf'
include { FEATURECOUNTS as FEATURECOUNTS_PROMOTER } from '../modules/featurecounts.nf'
include { FEATURECOUNTS as FEATURECOUNTS_GENEBODY } from '../modules/featurecounts.nf'
include { FEATURECOUNTS_MERGE as MERGE_PROMOTER   } from '../modules/featurecounts_merge.nf'
include { FEATURECOUNTS_MERGE as MERGE_GENEBODY   } from '../modules/featurecounts_merge.nf'

workflow quantification {
    take:
    bam_ch            // channel: tuple val(meta), path(bam) — one BAM per sample
    promoter_saf           // channel: path(promoter.saf)
    genebody_saf      // channel: path(genebody.saf)

    main:
    // PRO-seq 定量只取 read1（R2 是 5' 接头侧、无 Pol II 信号）：抽 R1-only
    // 单端 BAM 后再喂 featureCounts，避免 -p 把 R2 也计入。
    r1_ch = EXTRACT_R1(bam_ch)

    FEATURECOUNTS_PROMOTER(r1_ch, promoter_saf, 'promoter')
    FEATURECOUNTS_GENEBODY(r1_ch, genebody_saf, 'genebody')

    // Collect each region's per-sample count files -> merge into one matrix.
    MERGE_PROMOTER(FEATURECOUNTS_PROMOTER.out.counts.collect(), 'promoter')
    MERGE_GENEBODY(FEATURECOUNTS_GENEBODY.out.counts.collect(), 'genebody')

    emit:
    promoter_counts = FEATURECOUNTS_PROMOTER.out.counts     // per-sample raw counts
    genebody_counts = FEATURECOUNTS_GENEBODY.out.counts     // per-sample raw counts
    promoter_matrix = MERGE_PROMOTER.out.matrix             // merged matrix (gene_id, length, samples)
    genebody_matrix = MERGE_GENEBODY.out.matrix
}
