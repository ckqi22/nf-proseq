#!/usr/bin/env nextflow
//
// SUBWORKFLOW: prepare_genome
// Chains: Index build -> get chrom size -> gtf2bed
// Purpose: 
//

include { } from '../modules/fastqc.nf'
include { } from '../modules/cutadapt.nf'
include { } from '../modules/process_fastp.nf'

workflow prepare_genome {
    take:


    main:
    // ------------------------------------------------------------------
    // Step 1: (fastqc wrapper)
    // ------------------------------------------------------------------

    emit:

}
