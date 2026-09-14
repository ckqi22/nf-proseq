#!/usr/bin/env nextflow
//
// SUBWORKFLOW: align_bowtie2
// Chains: 
// Purpose: 
//

include { ALIGN as BOWTIE2_ALIGN} from '../modules/bowtie2/align.nf'
include { STAT                  } from '../modules/samtools/stat.nf'

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
    STAT(bam_bai_ch, fasta_ch)

    emit:
    bam         = BOWTIE2_ALIGN.out.bam
    bai         = BOWTIE2_ALIGN.out.bai
    bam_bai     = bam_bai_ch                        // tuple(meta, bam, bai) — 供 SPIKEIN_COUNT(idxstats 需 .bai)
    alignRate   = BOWTIE2_ALIGN.out.alignRate
    flagstat    = STAT.out.flagstat
    idxstats    = STAT.out.idxstats
    stats       = STAT.out.stats
}
