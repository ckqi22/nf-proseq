#!/usr/bin/env nextflow

include { GOKEGG } from '../modules/enrich/gokegg.nf'
include { GSEA   } from '../modules/enrich/gsea.nf'

workflow enrich {
    take:
    diff_dir

    main:
    GOKEGG(diff_dir)

    GSEA(diff_dir)

    emit:
    gokegg_result   = GOKEGG.out.result
    gsea_result     = GSEA.out.result
}
