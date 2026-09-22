#!/usr/bin/env nextflow
//
// SUBWORKFLOW: spikein
// spike-in 计数 + 缩放因子：
//   SPIKEIN_COUNT  从「主 + spike 合并参考」比对的 BAM 里按 spike 染色体统计每样本 read1 数
//                  （samtools view -f 64 抽 read1，SE 全量；无需 .bai）
//   SPIKEIN_SCALE  汇总全样本 count → factor = 1e6/spike_count 的因子表（bin/spikein_scale.R）
//
// 因子表下游用途：
//   - bigWig：pol2_count.nf 用 meta.spike_count 算 scale_spike → genomecov 出 _spike bigWig（与 _cpm 并列，不覆盖）
//   - 矩阵：bin/normalize.R --spike_factors 追加 .Spike 列（= count × factor），与 .Cpm/.Fpkm/.Rpkm 并列（不换分母）
//
//   占比质控：spike_count / total_mapped <= params.spike_min_fraction 时在 main: 直接 fail（默认 0.0 = 仅拒绝 0 spike reads）。
//

include { SPIKEIN_COUNT } from '../modules/spikein/spikein_count.nf'
include { SPIKEIN_SCALE } from '../modules/spikein/spikein_scale.nf'

workflow spikein {
    take:
    bam_bai       // channel: tuple(meta, bam, bai) — 合并参考比对的 BAM + 索引
    spike_chroms  // channel: path(spike.chroms) — SPIKEIN_CONCAT 产出
    total_mapped  // channel: tuple(meta, sample.total_mapped.txt) — read1 mapped 计数，用于 spike 占比质控

    main:
    def min_fraction = (params.spike_min_fraction ?: 0.0) as double

    counts = SPIKEIN_COUNT(bam_bai.map { meta, bam, _bai -> [meta, bam] }, spike_chroms).counts
        .join(total_mapped, by: [0])
        .map { meta, spike_file, mapped_file ->
            def spike_count         = spike_file.text.trim().toInteger()
            def total_mapped_count  = mapped_file.text.trim().toInteger()
            def spike_fraction      = spike_count.toDouble() / total_mapped_count
            if (spike_fraction <= min_fraction) {
                error "sample '${meta.sample}' spike fraction ${spike_fraction} (spike_count=${spike_count} / total_mapped=${total_mapped_count}) <= spike_min_fraction ${min_fraction}"
            }
            [meta, spike_file]
        }

    factors = SPIKEIN_SCALE(
        counts.map { meta, f -> [meta.sample, f] }        // 每样本 → [样本名, 文件]
              .toSortedList()                             // Nextflow 原生算子：收集成 list 并按样本名自然排序
              .map { list -> list.collect { it[1] } }     // 只留文件，顺序已固定
    ).factors

    emit:
    counts  = counts    // tuple(meta, spike_count.txt) — 每样本单行整数
    factors = factors   // path(spikein_scale_factors.tsv)
}
