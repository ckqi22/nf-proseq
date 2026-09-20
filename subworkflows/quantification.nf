#!/usr/bin/env nextflow
//
// SUBWORKFLOW: quantification
// Per-sample featureCounts over the gene body union annotation (one region type),
// then merge per-sample count files into one clean count matrix.
//
// featureCounts runs once per sample — the -p flag is decided per sample, so
// single-end and paired-end samples coexist correctly in one run.
//
// read1 单端化已上移到 align_bowtie2（EXTRACT_R1 每样本只抽一次），
// 此处直接吃现成的 r1_bam，不再自行抽 read1。
//

include { FEATURECOUNTS      } from '../modules/featurecounts.nf'
include { FEATURECOUNTS_MERGE } from '../modules/featurecounts_merge.nf'

workflow quantification {
    take:
    r1_bam             // channel: tuple val(meta), path(r1.bam) — read1 单端化 BAM（align_bowtie2 产出）
    genebody_union_saf // channel: tuple val(name), path(genebody_union.saf) — 区域名随文件走

    main:
    // PRO-seq 定量只取 read1（R2 是 5' 接头侧、无 Pol II 信号）：R1 已在 align_bowtie2
    // 抽好并剥 paired flag，直接喂 featureCounts，避免 -p 把 R2 也计入。
    FEATURECOUNTS(r1_bam, genebody_union_saf)

    // Collect per-sample count files -> merge into one matrix。
    FEATURECOUNTS_MERGE(FEATURECOUNTS.out.counts.collect(), 'genebody')

    emit:
    genebody_counts = FEATURECOUNTS.out.counts     // per-sample raw counts
    genebody_matrix = FEATURECOUNTS_MERGE.out.matrix   // merged matrix (gene_id, length, samples)
}
