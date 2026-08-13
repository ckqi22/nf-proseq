#!/usr/bin/env nextflow
//
// SUBWORKFLOW: align_bowtie2
// Chains: 
// Purpose: 
//

include { BOWTIE2_ALIGN } from '../modules/bowtie2/align.nf'
include { SAMTOOLS_STAT } from '../modules/samtools_stat.nf'

workflow align_bowtie2{
    take:
    reads_ch    // channel: [ val(meta), [ reads ] ]
    index_ch    // 
    fasta_ch    //

    main:
    // ------------------------------------------------------------------
    // Step 1: alignment
    // ------------------------------------------------------------------
    BOWTIE2_ALIGN(reads_ch, index_ch)

    bam_bai_ch = BOWTIE2_ALIGN.out.bam
        .join(BOWTIE2_ALIGN.out.bai, by: [0], remainder: true)
        .map {
            meta, bam, bai -> [meta, bam, bai]
        }    

    // ------------------------------------------------------------------
    // Step 2: samtools stat
    // ------------------------------------------------------------------    
    SAMTOOLS_STAT(bam_bai_ch, fasta_ch)

    emit:
    bam             = BOWTIE2_ALIGN.out.bam
    bai             = BOWTIE2_ALIGN.out.bai
    genomeRate      = BOWTIE2_ALIGN.out.genomeRate
    alignment_log   = BOWTIE2_ALIGN.out.alignment_log
    flagstat        = SAMTOOLS_STAT.out.flagstat
    idxstats        = SAMTOOLS_STAT.out.idxstats
    stats           = SAMTOOLS_STAT.out.stats
}
