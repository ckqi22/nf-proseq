#!/usr/bin/env nextflow
//
// SUBWORKFLOW: prepare_genome
// Resolve reference genome (fasta + bowtie2 index) and generate
// BED (tss/promoter/genebody) + tx2gene map + genebody-union SAF annotations
// from the reference GTF.
//

include { GTF2SAF } from '../modules/gtf2saf.nf'
include { GTF2BED } from '../modules/gtf2bed.nf'
include { BOWTIE2_BUILD } from '../modules/bowtie2/build.nf'
include { SPIKEIN_CONCAT } from '../modules/spikein/spikein_concat.nf'

workflow prepare_genome {
    take:
    config_ch         // channel: val(map) — genome_fasta / bowtie2_index / gtf / build / spike_fasta

    main:
    // ------------------------------------------------------------------
    // annotation BED：各 protein_coding transcript 的 TSS + promoter/genebody 窗口 + tx2gene
    // ------------------------------------------------------------------
    GTF2BED(config_ch.map { it -> it.gtf })

    // ------------------------------------------------------------------
    // annotation SAF：genebody union（featureCounts 定量）
    // ------------------------------------------------------------------
    GTF2SAF(config_ch.map { it -> it.gtf })

    // ------------------------------------------------------------------
    // reference: fasta + bowtie2 index
    // ------------------------------------------------------------------
    // `index` 统一为「索引所在目录」（path，stage 后 align 从 *.1.bt2 反推前缀）。
    //   - 无 spike：直接用 DB 预建索引的父目录（bowtie2_index 是无后缀前缀）
    //   - 有 spike + spike_index：用预构建合并索引的父目录（跳过 concat+build）
    //   - 有 spike 无 spike_index：SPIKEIN_CONCAT 合并主+spike fasta → BOWTIE2_BUILD 现建合并索引
    // 注释（BED/SAF/GTF）只用主 GTF，spike 染色体不进 SAF/BED，不污染计数。
    fasta = config_ch.map { it -> [ [id: it.build], it.genome_fasta ] }

    // chrom.sizes（主基因组，无 spike）：从参考 FASTA 的 .fai 索引现取（samtools faidx 产物，
    // 染色体集 = 实际比对基因组，不含 DB genomesizefile 里的 non-m6a/m6a 等额外 contig）。
    // 5 列 .fai 由下游 BIGWIGAVERAGE cut -f1,2 转成 2 列 chrom.sizes。
    chrom_sizes = config_ch.map { it -> file("${it.genome_fasta}.fai") }

    def spike_split = config_ch.branch {
        has_spike: (it.spike_fasta ?: '').trim() || (it.spike_index ?: '').trim()
        no_spike:  !((it.spike_fasta ?: '').trim() || (it.spike_index ?: '').trim())
    }

    index_prebuilt = spike_split.no_spike.map { it -> [ [id: it.build], file(it.bowtie2_index).parent ] }

    // has_spike 再分：给了 spike_index → 用预构建合并索引（跳过 concat+build）；否则 concat+build 现建
    def spike_mode = spike_split.has_spike.branch {
        prebuilt: (it.spike_index ?: '').trim()
        build:    !((it.spike_index ?: '').trim())
    }
    index_spike_prebuilt = spike_mode.prebuilt.map { it -> [ [id: it.build], file(it.spike_index).parent ] }

    concat_ch = SPIKEIN_CONCAT(
        spike_mode.build.map { it -> it.genome_fasta },
        spike_mode.build.map { it -> it.spike_fasta })
    built = BOWTIE2_BUILD(concat_ch.combined_fasta.map { f -> [ [id: 'combined'], f ] }).index

    index = index_prebuilt.mix(index_spike_prebuilt, built).first()

    // spike 染色体名单：预构建时用现成文件；现建时用 CONCAT 产出。加 .first() 保持广播（同 index）
    spike_chroms = spike_mode.prebuilt.map { it -> file(it.spike_chroms) }
                     .mix(concat_ch.spike_chroms).first()

    emit:
    index               = index
    spike_chroms        = spike_chroms
    fasta               = fasta
    chrom_sizes         = chrom_sizes
    tx2gene             = GTF2BED.out.tx2gene               // transcript_id → gene_id 映射 → merge 加 gene_id 列
    tss_bed             = GTF2BED.out.tss_bed               // transcript TSS BED6 → TSS metagene
    promoter_bed        = GTF2BED.out.promoter_bed          // transcript promoter → pol2_count 单碱基 promoter 计数
    genebody_bed        = GTF2BED.out.genebody_bed          // transcript genebody → pol2_count 单碱基 genebody 计数
    genebody_union_saf  = GTF2SAF.out.genebody_union_saf    // 所有 transcript genebody union → quantification(featureCounts)
}
