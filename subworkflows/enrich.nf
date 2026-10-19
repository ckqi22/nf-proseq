#!/usr/bin/env nextflow
//
// SUBWORKFLOW: enrich
// 差异结果的 GO/KEGG 富集与 GSEA：
//   GOKEGG  enrichment.R --methods NORMAL（GO + KEGG 富集）
//   GSEA    enrichment.R --methods GSEA_KEGG（GSEA）
// 输入为通过 CHECKDE 的差异结果目录（diff.out.passed）；两者共享同一diff_dir。
//

include { GOKEGG } from '../modules/enrich/gokegg.nf'
include { GSEA   } from '../modules/enrich/gsea.nf'

workflow enrich {
    take:
    diff_dir  // channel: path — 通过 CHECKDE 的差异结果目录（diff.out.passed）

    main:
    GOKEGG(diff_dir)

    GSEA(diff_dir)

    emit:
    gokegg_result   = GOKEGG.out.result   // path(gokegg_result/*) — GO/KEGG 富集结果
    gsea_result     = GSEA.out.result     // path(gsea_result/*) — GSEA_KEGG 结果
}
