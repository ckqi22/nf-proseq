#!/usr/bin/env nextflow
//
// SUBWORKFLOW: preprocess
// Chains: fastc -> fastp -> process_fastp
// Purpose: Read QC, adapter trimming
//

include { FASTQC as FASTQC_RAW     } from '../modules/fastqc.nf'
include { FASTP                    } from '../modules/fastp.nf'
include { FASTQC as FASTQC_TRIMMED } from '../modules/fastqc.nf'
include { PROCESS_FASTP            } from '../modules/process_fastp.nf'

workflow preprocess {
    take:
    read_ch                          // channel: tuple val(meta), path(reads)
    adapter_type                     // val: adapter type (I, UMI, HT, SP)

    main:
    // ------------------------------------------------------------------
    // Step 1: (fastqc wrapper)
    // ------------------------------------------------------------------
    FASTQC_RAW(read_ch)

    // ------------------------------------------------------------------
    // Step 2: Adapter trimming and read QC with cutadapt (fastp wrapper)
    // ------------------------------------------------------------------
    FASTP(read_ch, adapter_type)

    // ------------------------------------------------------------------
    // Step 3: (fastqc wrapper)
    // ------------------------------------------------------------------
    FASTQC_TRIMMED(FASTP.out.trimmed_reads)

    // ------------------------------------------------------------------
    // Step 4: Process fastp output (base quality plots, read statistics)
    // ------------------------------------------------------------------
    PROCESS_FASTP(FASTP.out.json.map { _meta, json -> json }.collect())

    emit:
    // Cutadapt / QC outputs
    fastqc_raw_zip      = FASTQC_RAW.out.zip
    fastqc_raw_html     = FASTQC_RAW.out.html
    fastp_json          = FASTP.out.json
    fastp_html          = FASTP.out.html
    fastqc_trimmed_zip  = FASTQC_TRIMMED.out.zip
    fastqc_trimmed_html = FASTQC_TRIMMED.out.html    
    fastp_log           = FASTP.out.log
    trimmed_reads       = FASTP.out.trimmed_reads
    base_quality_plot   = PROCESS_FASTP.out.base_quality_plot
    statistics          = PROCESS_FASTP.out.statistics
}
