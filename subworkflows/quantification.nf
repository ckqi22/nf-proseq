#!/usr/bin/env nextflow
//
// SUBWORKFLOW: quantification
// Two regions per gene: TSS pause / gene body
// featureCounts runs once per region with all samples' BAMs together.
//

include { FEATURECOUNTS as FEATURECOUNTS_TSS } from '../modules/featurecounts.nf'
include { FEATURECOUNTS as FEATURECOUNTS_GENEBODY } from '../modules/featurecounts.nf'

workflow quantification {
    take:
    bam_ch            // channel: tuple val(meta), path(bams) — all samples' BAMs (collected in main.nf)
    tss_saf           // channel: path(tss.saf)
    genebody_saf      // channel: path(genebody.saf)

    main:
    FEATURECOUNTS_TSS(bam_ch, tss_saf, 'tss')
    FEATURECOUNTS_GENEBODY(bam_ch, genebody_saf, 'genebody')

    emit:
    tss_counts      = FEATURECOUNTS_TSS.out.counts
    genebody_counts = FEATURECOUNTS_GENEBODY.out.counts
}
