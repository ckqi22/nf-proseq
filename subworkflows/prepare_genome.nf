#!/usr/bin/env nextflow
//
// SUBWORKFLOW: prepare_genome
// Resolve reference genome (fasta + bowtie2 index) and generate
// TSS + gene body SAF annotations from the reference GTF.
//

include { GTF2SAF } from '../modules/gtf2saf.nf'

workflow prepare_genome {
    take:
    config_ch         // channel: val(map) — genome_fasta / bowtie2_index / gtf / build

    main:
    // ------------------------------------------------------------------
    // annotation: TSS + gene body SAF
    // ------------------------------------------------------------------
    GTF2SAF(config_ch.map { it -> it.gtf })

    // ------------------------------------------------------------------
    // reference: fasta (single file) + bowtie2 index (prebuilt prefix)
    // ------------------------------------------------------------------
    // `index` is the extensionless bowtie2 prefix (absolute path in the DB);
    // passed as a val so align runs `bowtie2 -x ${index}` without staging.
    index = config_ch.map { it -> [ [id: it.build], it.bowtie2_index ] }
    fasta = config_ch.map { it -> [ [id: it.build], it.genome_fasta ] }

    // TODO(on-demand build): when bowtie2_index is empty (new species), build
    //   from genome_fasta via BOWTIE2_BUILD (modules/bowtie2/build.nf).
    // TODO(spike-in): concat genome + spike fasta/gtf into a single reference
    //   before SAF/index generation.

    emit:
    index           = index
    fasta           = fasta
    tss_saf         = GTF2SAF.out.tss_saf
    genebody_saf    = GTF2SAF.out.genebody_saf
}
