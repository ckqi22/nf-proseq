#!/usr/bin/env nextflow
//
// SUBWORKFLOW: spikein
// spike-in 计数 + 缩放因子：
//   SPIKEIN_COUNT  从「主 + spike 合并参考」比对的 BAM 里按 spike 染色体统计每样本 read 数
//                  （samtools idxstats，需 .bai）
//   SPIKEIN_SCALE  汇总全样本 count → factor = 1e6/spike_count 的因子表（bin/spikein_scale.R）
//
// 因子表下游用途：
//   - bigWig：genomecov.nf 就地用 meta.spike_count 缩放（见 modules/bedtools/genomecov.nf）
//   - 矩阵：bin/normalize.R --spike_factors 把 CPM/FPKM 分母换成 spike_count
//

include { SPIKEIN_COUNT } from '../modules/spikein/spikein_count.nf'
include { SPIKEIN_SCALE } from '../modules/spikein/spikein_scale.nf'

workflow spikein {
    take:
    bam_bai       // channel: tuple(meta, bam, bai) — 合并参考比对的 BAM + 索引
    spike_chroms  // channel: path(spike.chroms) — SPIKEIN_CONCAT 产出（spike fasta 头第一 token）

    main:
    counts = SPIKEIN_COUNT(bam_bai, spike_chroms).counts

    factors = SPIKEIN_SCALE(counts.map { _m, f -> f }.collect()).factors

    emit:
    counts  = counts    // tuple(meta, spike_count.txt) — 每样本单行整数
    factors = factors   // path(spikein_scale_factors.tsv)
}
