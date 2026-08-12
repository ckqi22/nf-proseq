#!/usr/bin/env nextflow
//
// SUBWORKFLOW: metagene_analysis
// Chains: metagene_tss + metagene_tes (strand-specific metagene profiles)
//

include { metagene_tss } from '../modules/metagene_tss.nf'
include { metagene_tes } from '../modules/metagene_tes.nf'

workflow metagene_analysis {
    take:
    plus_bam_list     // channel: path — file listing plus-strand BAM paths
    minus_bam_list    // channel: path — file listing minus-strand BAM paths
    config_ch         // channel: val(config)

    main:
    metagene_tss(plus_bam_list, minus_bam_list, config_ch)
    metagene_tes(plus_bam_list, minus_bam_list, config_ch)

    emit:
    tss_plus_pdf  = metagene_tss.out.tss_plus_pdf
    tss_minus_pdf = metagene_tss.out.tss_minus_pdf
    tes_plus_pdf  = metagene_tes.out.tes_plus_pdf
    tes_minus_pdf = metagene_tes.out.tes_minus_pdf
}
