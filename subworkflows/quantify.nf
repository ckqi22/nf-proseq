#!/usr/bin/env nextflow
//
// SUBWORKFLOW: quantify
// Chains: featurecounts_gb + featurecounts_tss (per-sample strand-specific)
// Purpose: PRO-seq gene body and TSS window quantification
//

include { featurecounts_gb  } from '../modules/featurecounts_gb.nf'
include { featurecounts_tss } from '../modules/featurecounts_tss.nf'

workflow quantify {
    take:
    strand_bams       // channel: tuple val(meta), path(plus_bam), path(minus_bam)
    config_ch          // channel: val(config)

    main:
    // Build input channel: tuple val(meta), path(plus_bam), path(minus_bam), val(config)
    quant_input_ch = strand_bams
        .combine(config_ch)
        .map { meta, plus_bam, minus_bam, cfg ->
            [meta, plus_bam, minus_bam, cfg]
        }

    // Gene body quantification (TSS+301bp to TES, exon-level, strand-specific)
    featurecounts_gb(quant_input_ch)

    // TSS window quantification (-50bp to +300bp, strand-specific)
    featurecounts_tss(quant_input_ch)

    emit:
    gene_body_plus   = featurecounts_gb.out.plus_counts
    gene_body_minus  = featurecounts_gb.out.minus_counts
    gene_body_counts = featurecounts_gb.out.combined_counts

    tss_plus_counts  = featurecounts_tss.out.tss_plus_counts
    tss_minus_counts = featurecounts_tss.out.tss_minus_counts
    tss_counts       = featurecounts_tss.out.combined_tss_counts
}
