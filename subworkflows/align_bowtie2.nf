#!/usr/bin/env nextflow
//
// SUBWORKFLOW: align_bowtie2
// Chains: alignment -> read1 单端化 -> samtools stat
// Purpose: 产出全量 bam(+bai)、read1 单端 bam(r1_bam)、比对统计
//

include { ALIGN as BOWTIE2_ALIGN} from '../modules/bowtie2/align.nf'
include { STAT                  } from '../modules/samtools/stat.nf'
include { EXTRACT_R1            } from '../modules/samtools/extract_r1.nf'

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

    // ------------------------------------------------------------------
    // Step 2: read1 单端化（PRO-seq 只取 read1；SE 直接 cp 原 bam）
    //   r1_bam 供下游 quantification(featureCounts) 与 pol2_count(GENOMECOV) 共享，
    //   每样本只抽一次，避免两处各自 inline 抽 read1。
    // ------------------------------------------------------------------
    EXTRACT_R1(BOWTIE2_ALIGN.out.bam)

    bam_bai_ch = BOWTIE2_ALIGN.out.bam
        .join(BOWTIE2_ALIGN.out.bai, by: [0], remainder: true)
        .map {
            meta, bam, bai -> [meta, bam, bai]
        }

    // ------------------------------------------------------------------
    // Step 3: samtools stat
    // ------------------------------------------------------------------
    STAT(bam_bai_ch, fasta_ch)

    emit:
    bam         = BOWTIE2_ALIGN.out.bam
    bai         = BOWTIE2_ALIGN.out.bai
    bam_bai     = bam_bai_ch                        // tuple(meta, bam, bai) — 供 SPIKEIN_COUNT(idxstats 需 .bai)
    r1_bam      = EXTRACT_R1.out.r1_bam             // tuple(meta, r1.bam) — read1 单端化 BAM（SE=原 bam 拷贝）
    total_mapped = EXTRACT_R1.out.total_mapped       // tuple(meta, sample.total_mapped.txt) — read1 mapped count（归一化分母）
    alignRate   = BOWTIE2_ALIGN.out.alignRate
    flagstat    = STAT.out.flagstat
    idxstats    = STAT.out.idxstats
    stats       = STAT.out.stats
}
