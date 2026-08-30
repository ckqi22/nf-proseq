#!/usr/bin/env nextflow
//
// SUBWORKFLOW: quantification
// Per-sample featureCounts over the gene body union annotation (one region type),
// then merge per-sample count files into one clean count matrix.
//
// featureCounts runs once per sample — the -p flag is decided per sample, so
// single-end and paired-end samples coexist correctly in one run.
//

include { EXTRACT_R1         } from '../modules/samtools/extract_r1.nf'
include { FEATURECOUNTS      } from '../modules/featurecounts.nf'
include { FEATURECOUNTS_MERGE } from '../modules/featurecounts_merge.nf'

workflow quantification {
    take:
    bam_ch             // channel: tuple val(meta), path(bam) — one BAM per sample
    genebody_union_saf // channel: tuple val(name), path(genebody_union.saf) — 区域名随文件走

    main:
    // PRO-seq 定量只取 read1（R2 是 5' 接头侧、无 Pol II 信号）：抽 R1-only
    // 单端 BAM 后再喂 featureCounts，避免 -p 把 R2 也计入。
    r1_ch = EXTRACT_R1(bam_ch)

    FEATURECOUNTS(r1_ch, genebody_union_saf)

    // Collect per-sample count files -> merge into one matrix。
    FEATURECOUNTS_MERGE(FEATURECOUNTS.out.counts.collect(), 'genebody')

    emit:
    genebody_counts = FEATURECOUNTS.out.counts     // per-sample raw counts
    genebody_matrix = FEATURECOUNTS_MERGE.out.matrix   // merged matrix (gene_id, length, samples)
}
