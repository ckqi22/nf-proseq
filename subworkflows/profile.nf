#!/usr/bin/env nextflow
//
// SUBWORKFLOW: profile
// 把一张 count 矩阵（gene_id, length, samples）归一化（CPM/FPKM）并拼上本地基因注释，
// 输出一张 profile 表：Gene_id  .Count  .Fpkm  .Cpm  <annotation cols...>。
// 对 featureCounts（promoter / genebody）与 Pol II 单碱基矩阵通用。
//

include { NORMALIZE } from '../modules/normalize.nf'
include { ANNOTATE  } from '../modules/annotate.nf'

workflow profile {
    take:
    matrix        // channel: path(raw count matrix)
    annotation    // channel: val(string) — gene_annotation 路径（同 diff 的 annotation_ch）
    methods       // channel: val ("cpm,fpkm")
    name          // channel: val ("genebody" | "promoter" | "pol2")

    main:
    NORMALIZE(matrix, methods, name)
    ANNOTATE(matrix, NORMALIZE.out.normalized, annotation, name)

    emit:
    annotated = ANNOTATE.out.annotated   // 每矩阵一张表
}
