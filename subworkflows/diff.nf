#!/usr/bin/env nextflow
//
// SUBWORKFLOW: diff
// Differential analysis on gene body counts.
//   Method-agnostic take/emit contract. Currently only DESeq2 is wired; the
//   dispatch point below is where an edgeR (or other) method can be added.
//

include { DESEQ2 } from '../modules/diff/deseq2.nf'
// TODO(edgeR): include { EDGER } from '../modules/diff/edger.nf'

workflow diff {
    take:
    genebody_matrix    // channel: path — gene body count matrix (gene_id, length, samples)
    groups             // channel: val(map) — group_name -> [samples] (from samplesheet)
    annotation         // channel: path — gene annotation (gene_id, gene_name, gene_biotype)

    main:
    // Build the differential config YAML (group + compared-groups). Shared by
    // DESeq2 and a future edgeR method, so it lives here in the dispatch layer.
    config_yml = groups.map { g ->
        def lines = ["group:"]
        g.each { name, samples -> lines << "  ${name}: \"${samples.join(', ')}\"" }
        lines << "compared-groups:"
        (params.compared_groups ?: []).each { e -> lines << "  - \"${e}\"" }
        lines.join('\n')
    }

    // Method dispatch — currently DESeq2 only.
    // TODO(edgeR): route to EDGER (or others) here, e.g. on params.diff.tool.
    DESEQ2(genebody_matrix, config_yml, annotation)

    emit:
    results = DESEQ2.out.results
}
