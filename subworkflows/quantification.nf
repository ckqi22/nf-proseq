#!/usr/bin/env nextflow
//
// SUBWORKFLOW: quantification
// Three regions per gene: full gene / TSS pause / gene body
// All in one output file per sample.
//

include {  }
include { FEATURECOUNTS } from '../modules/featurecounts.nf'

workflow quantification {
    take:
    bam_ch            // channel: tuple val(meta), path(bam)
    config_ch         // channel: val(config)

    main:
    quant_input_ch = bam_ch
        .combine(config_ch)
        .map { meta, bam, cfg -> [meta, bam, cfg] }

    featurecounts_proseq(quant_input_ch)

    emit:
    counts  = featurecounts_proseq.out.counts
    summary = featurecounts_proseq.out.summary
}
