#!/usr/bin/env nextflow
//
// SUBWORKFLOW: metagene
// TSS-centered metagene analysis from strand-specific CPM bigWigs:
//   - PLOT_PROFILE: 平均 TSS metagene 曲线（所有基因取均值）
//   - PLOT_HEATMAP: 逐基因 TSS 热图（每行一个基因，按窗口内信号排序）
// 两者共用同一份 tss.bed（代表 transcript TSS BED6，由 prepare_genome 生成）。
//

include { PLOT_PROFILE } from '../modules/metagene/plot_profile.nf'
include { PLOT_HEATMAP } from '../modules/metagene/plot_heatmap.nf'

workflow metagene {
    take:
    bigwig_cpm   // channel: tuple(meta, plus_bw, minus_bw) — CPM bigWigs
    tss_bed      // channel: path(tss.bed)

    main:
    PLOT_PROFILE(bigwig_cpm, tss_bed)
    PLOT_HEATMAP(bigwig_cpm, tss_bed)

    emit:
    plot   = PLOT_PROFILE.out.plot.mix(PLOT_HEATMAP.out.plot)
    matrix = PLOT_PROFILE.out.matrix.mix(PLOT_HEATMAP.out.matrix)
}
