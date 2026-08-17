#!/usr/bin/env nextflow
//
// SUBWORKFLOW: quantification
// Per-sample featureCounts over two regions per gene (TSS pause / gene body),
// then merge each region's per-sample count files into one clean count matrix.
//
// featureCounts runs once per sample per region — the -p flag is decided per
// sample, so single-end and paired-end samples coexist correctly in one run.
//

include { FEATURECOUNTS as FEATURECOUNTS_TSS      } from '../modules/featurecounts.nf'
include { FEATURECOUNTS as FEATURECOUNTS_GENEBODY } from '../modules/featurecounts.nf'
include { FEATURECOUNTS_MERGE as MERGE_TSS        } from '../modules/featurecounts_merge.nf'
include { FEATURECOUNTS_MERGE as MERGE_GENEBODY   } from '../modules/featurecounts_merge.nf'

workflow quantification {
    take:
    bam_ch            // channel: tuple val(meta), path(bam) — one BAM per sample
    tss_saf           // channel: path(tss.saf)
    genebody_saf      // channel: path(genebody.saf)

    main:
    FEATURECOUNTS_TSS(bam_ch, tss_saf, 'tss')
    FEATURECOUNTS_GENEBODY(bam_ch, genebody_saf, 'genebody')

    // Collect each region's per-sample count files -> merge into one matrix.
    MERGE_TSS(FEATURECOUNTS_TSS.out.counts.collect(), 'tss')
    MERGE_GENEBODY(FEATURECOUNTS_GENEBODY.out.counts.collect(), 'genebody')

    emit:
    tss_counts       = FEATURECOUNTS_TSS.out.counts        // per-sample raw counts
    genebody_counts  = FEATURECOUNTS_GENEBODY.out.counts   // per-sample raw counts
    tss_matrix       = MERGE_TSS.out.matrix                // merged matrix (gene_id, length, samples)
    genebody_matrix  = MERGE_GENEBODY.out.matrix
}
