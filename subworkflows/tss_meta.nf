#!/usr/bin/env nextflow
//
// SUBWORKFLOW: tss_meta
// TSS metagene profile via deepTools:
//   bigWig (来自 pol2_profile, gene-strand) -> computeMatrix (reference-point TSS) -> plotProfile
//   正链基因、负链基因分别出图：各自一张 TSS profile PDF。
//

include { COMPUTEMATRIX } from '../modules/deeptools/computeMatrix.nf'
include { PLOTPROFILE   } from '../modules/deeptools/plotProfile.nf'

workflow tss_meta {
    take:
    bigwig          // channel: tuple(meta, plus_bw, minus_bw)  gene-strand：plus=+基因信号、minus=-基因信号（取负）
    gene_bed

    main:
    COMPUTEMATRIX(bigwig, gene_bed)

    PLOTPROFILE(COMPUTEMATRIX.out.matrix)

    emit:
    matrix  = COMPUTEMATRIX.out.matrix          // tuple(meta, plus_matrix, minus_matrix)
    profile = PLOTPROFILE.out.profile           // 两张 PDF：<sample>_plus/minus_TSS_meta.pdf
}
