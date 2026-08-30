#!/usr/bin/env nextflow
//
// SUBWORKFLOW: prepare_genome
// Resolve reference genome (fasta + bowtie2 index) and generate
// promoter + genebody SAF annotations from the reference GTF.
//

include { GTF2SAF } from '../modules/gtf2saf.nf'

workflow prepare_genome {
    take:
    config_ch         // channel: val(map) — genome_fasta / bowtie2_index / gtf / build

    main:
    // ------------------------------------------------------------------
    // annotation: promoter + genebody SAF
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
    index               = index
    fasta               = fasta
    promoter_bed        = GTF2SAF.out.promoter_bed        // 最长 transcript promoter → pol2_count 单碱基 promoter 计数
    genebody_bed        = GTF2SAF.out.genebody_bed        // 最长 transcript genebody → pol2_count 单碱基 genebody 计数
    genebody_union_saf  = GTF2SAF.out.genebody_union_saf  // 所有 transcript genebody union → quantification(featureCounts)
    gene_bed            = GTF2SAF.out.gene_bed            // gene 级跨度 → tss_meta + SIGNAL_TABLE(信号表 intersect)
    // plus_genes_bed  = GTF2SAF.out.plus_genes_bed    // 未使用，先注释掉
    // minus_genes_bed = GTF2SAF.out.minus_genes_bed   // 未使用，先注释掉
}
