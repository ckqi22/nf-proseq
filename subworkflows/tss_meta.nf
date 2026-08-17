#!/usr/bin/env nextflow
//
// SUBWORKFLOW: tss_meta
// TSS metagene profile via deepTools:
//   BAM -> bigWig (bamCoverage) -> computeMatrix (reference-point TSS) -> plotProfile
//

include { BAMCOVERAGE   } from '../modules/deeptools/bamcoverage.nf'
include { COMPUTEMATRIX } from '../modules/deeptools/computeMatrix.nf'
include { PLOTPROFILE   } from '../modules/deeptools/plotProfile.nf'

workflow tss_meta {
    take:
    bam    // channel: tuple(meta, bam)
    bai    // channel: tuple(meta, bai)
    gtf    // channel: val(gtf_path)

    main:
    BAMCOVERAGE(bam.join(bai))

    COMPUTEMATRIX(BAMCOVERAGE.out.bigwig, gtf)

    PLOTPROFILE(COMPUTEMATRIX.out.matrix)

    emit:
    bigwig  = BAMCOVERAGE.out.bigwig
    matrix  = COMPUTEMATRIX.out.matrix
    profile = PLOTPROFILE.out.profile
}
