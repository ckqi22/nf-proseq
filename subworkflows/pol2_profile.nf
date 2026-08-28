#!/usr/bin/env nextflow
//
// SUBWORKFLOW: pol2_profile
// Pol II 单碱基 5' 端活性位点覆盖(bedtools genomecov，gene-strand 约定):
//   GENOMECOV 从 BAM 生成单碱基 bedGraph + bigWig：
//     _plus = + 链基因信号（正）；_minus = - 链基因信号（取负）
//   随后按基因链 sum 得到每基因 count 矩阵（+ 链基因读 _plus、- 链基因读 _minus 取负回正）。
//

include { GENOMECOV        } from '../modules/bedtools/genomecov.nf'
include { POL2_COUNT       } from '../modules/pol2_count.nf'
include { POL2_COUNT_MERGE } from '../modules/pol2_count_merge.nf'

workflow pol2_profile {
    take:
    bam         // channel: tuple(meta, bam)
    gene_bed    // channel: path(gene_bed)

    main:
    GENOMECOV(bam)

    POL2_COUNT(GENOMECOV.out.bedgraph, gene_bed)

    POL2_COUNT_MERGE(POL2_COUNT.out.counts.collect())

    emit:
    bigwig   = GENOMECOV.out.bigwig         // tuple(meta, plus_bw, minus_bw)  gene-strand
    bedGraph = GENOMECOV.out.bedgraph       // tuple(meta, plus_bg, minus_bg)
    counts   = POL2_COUNT.out.counts        // per-sample gene_id length count
    matrix   = POL2_COUNT_MERGE.out.matrix  // merged matrix (gene_id, length, samples)
}
