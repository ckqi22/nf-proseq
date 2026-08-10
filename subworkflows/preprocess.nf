#!/usr/bin/env nextflow
//
// SUBWORKFLOW: preprocess
// Chains: cutadapt -> process_cutadapt_output -> alignment_bowtie2 -> strand_split
// Purpose: Read QC, adapter trimming, alignment, and strand separation for PRO-seq data
// NOTE: config_ch is parsed by the caller (main.nf) and passed in — no duplicate parse_config call.
//

include { cutadapt                } from '../modules/cutadapt.nf'
include { process_cutadapt_output } from '../modules/process_cutadapt_output.nf'
include { alignment_bowtie2       } from '../modules/alignment_bowtie2.nf'
include { strand_split            } from '../modules/strand_split.nf'

workflow preprocess {
    take:
    read_ch                          // channel: tuple val(meta), path(reads)
    adapter_type                     // val: adapter type (I, UMI, HT, SP)
    config_ch                        // val: parsed config map (from parse_config in main.nf)

    main:
    // ------------------------------------------------------------------
    // Step 1: Adapter trimming and read QC with cutadapt (fastp wrapper)
    // ------------------------------------------------------------------
    cutadapt(read_ch, adapter_type)

    // ------------------------------------------------------------------
    // Step 2: Process cutadapt output (base quality plots, read statistics)
    // ------------------------------------------------------------------
    process_cutadapt_output(cutadapt.out.json)

    // ------------------------------------------------------------------
    // Step 3: Bowtie2 alignment (very-sensitive, --rf for dUTP strand-specific)
    // ------------------------------------------------------------------
    alignment_bowtie2(cutadapt.out.trimmed_reads, config_ch)

    // ------------------------------------------------------------------
    // Step 4: Strand separation
    // Forward reads (-F 0x10) go to plus.bam for plus-strand genes
    // Reverse reads (-f 0x10) go to minus.bam for minus-strand genes
    // ------------------------------------------------------------------
    strand_split(alignment_bowtie2.out.bam)

    emit:
    // Strand-split BAMs: tuple val(meta), path(plus_bam), path(minus_bam)
    strand_bams = strand_split.out.strand_bams

    // Unsplitt BAMs (for pol2_profiling and other processes needing combined BAM)
    bam    = alignment_bowtie2.out.bam
    bai    = alignment_bowtie2.out.bai
    genomeRate     = alignment_bowtie2.out.genomeRate
    alignment_log  = alignment_bowtie2.out.alignment_log

    // Strand indices
    plus_bai  = strand_split.out.plus_bai
    minus_bai = strand_split.out.minus_bai

    // Cutadapt / QC outputs
    cutadapt_json       = cutadapt.out.json
    cutadapt_html       = cutadapt.out.html
    cutadapt_log        = cutadapt.out.log
    trimmed_reads       = cutadapt.out.trimmed_reads
    base_quality_plot   = process_cutadapt_output.out.base_quality_plot
    statistics          = process_cutadapt_output.out.statistics
    statistics_log      = process_cutadapt_output.out.statistics_log
}
