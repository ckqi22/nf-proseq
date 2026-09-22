#!/usr/bin/env nextflow
//
// SUBWORKFLOW: diff
// Differential analysis on gene body counts.
//   Method-agnostic take/emit contract. Method selected by params.diff.tool
//   (deseq2 default / edger); downstream CHECKDE/plots/enrich are method-agnostic.
//

include { DESEQ2  } from '../modules/diff/deseq2.nf'
include { EDGER   } from '../modules/diff/edger.nf'
include { CHECKDE } from '../modules/diff/checkDE.nf'
include { DE_PLOT  } from '../modules/diff/de_plot.nf'
include { PCA_PLOT } from '../modules/diff/pca_plot.nf'

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
        (params.diff?.compared_groups ?: []).each { e -> lines << "  - \"${e}\"" }
        lines.join('\n')
    }

    // Method dispatch — params.diff.tool: deseq2 (default) | edger.
    def method = params.diff?.tool ?: 'deseq2'
    if (method == 'edger') {
        EDGER(genebody_matrix, config_yml, annotation)
        result_ch = EDGER.out.result
    } else {
        DESEQ2(genebody_matrix, config_yml, annotation)
        result_ch = DESEQ2.out.result
    }

    // DE count QC (soft-skip): CHECKDE always exits 0, writing a PASS/FAIL
    // flag into checkDE_result.txt. Route the diff_dir to `passed` only when
    // the check passes, so enrich runs exclusively on samples that passed.
    CHECKDE(result_ch)

    // DE visualization (runs on the raw result, independent of the CHECKDE
    // gate — plots are useful even when the DE count is below threshold).
    DE_PLOT(result_ch)
    PCA_PLOT(result_ch)

    def checkde_br = CHECKDE.out.checked.branch {
        pass: it[0].text.trim().startsWith("PASS")
        fail: !it[0].text.trim().startsWith("PASS")
    }

    emit:
    result         = result_ch                                        // unchanged — published to 06.
    passed         = checkde_br.pass.map { _flag, dir -> dir }        // for enrich (pass-branch only)
    checkde_result = CHECKDE.out.checked.map { flag, _dir -> flag }   // QC receipt — published to 06.
    de_plot        = DE_PLOT.out.plot
    pca_plot       = PCA_PLOT.out.pca
}
