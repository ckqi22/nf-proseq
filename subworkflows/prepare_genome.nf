#!/usr/bin/env nextflow
//
// SUBWORKFLOW: prepare_genome
// Resolve reference genome (fasta + bowtie2 index) and generate
// BED (tss/promoter/genebody/gene) + genebody-union SAF annotations
// from the reference GTF.
//

include { GTF2SAF } from '../modules/gtf2saf.nf'
include { LONGEST_TX } from '../modules/longest_tx.nf'
include { GTF2BED } from '../modules/gtf2bed.nf'

workflow prepare_genome {
    take:
    config_ch         // channel: val(map) — genome_fasta / bowtie2_index / gtf / build

    main:
    // ------------------------------------------------------------------
    // representative transcript (longest protein_coding per gene, exon-sum) → GTF
    // 单一来源：PI(分支A) 与 TSS metagene 均由此 GTF 派生。
    // ------------------------------------------------------------------
    LONGEST_TX(config_ch.map { it -> it.gtf })

    // ------------------------------------------------------------------
    // annotation BED：代表 transcript 的 TSS 碱基 + promoter/genebody 窗口 + gene 跨度
    // ------------------------------------------------------------------
    GTF2BED(config_ch.map { it -> it.gtf }, LONGEST_TX.out.representative_gtf)

    // ------------------------------------------------------------------
    // annotation SAF：genebody union（featureCounts 定量）
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
    representative_gtf  = LONGEST_TX.out.representative_gtf // 代表转录本 GTF（供 GTF2BED / 下游）
    tss_bed             = GTF2BED.out.tss_bed               // 代表 transcript TSS BED6 → TSS metagene
    promoter_bed        = GTF2BED.out.promoter_bed          // 代表 transcript promoter → pol2_count 单碱基 promoter 计数
    genebody_bed        = GTF2BED.out.genebody_bed          // 代表 transcript genebody → pol2_count 单碱基 genebody 计数
    genebody_union_saf  = GTF2SAF.out.genebody_union_saf    // 所有 transcript genebody union → quantification(featureCounts)
    gene_bed            = GTF2BED.out.gene_bed              // gene 级跨度 → SIGNAL_TABLE(信号表 intersect)
}
