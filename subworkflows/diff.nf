#!/usr/bin/env nextflow
//
// SUBWORKFLOW: diff
// 基因体计数的差异分析。
//   take/emit 契约与方法无关；方法由 params.diff.tool 选择（默认 edger / deseq2）；
//   下游 CHECKDE/图/enrich 均方法无关。
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

    // Method dispatch — params.diff.tool: edger (default) | deseq2.
    def method = params.diff?.tool ?: 'edger'
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

    // DE count QC（软跳过）：CHECKDE 的 flag 为 FAIL 时用 log.warn 提示（流程继续，仅跳过富集）。
    // 用最常用的 map/filter 实现：map 里算 pass 布尔并顺手 log.warn，filter 只保留 PASS 送下游。
    def checked = CHECKDE.out.checked.map { flag, dir ->
        def pass = flag.text.trim().startsWith("PASS")
        if (!pass) {
            flag.text.readLines()
                .findAll { line -> line.startsWith('Error:') }
                .each { line ->
                    log.warn "${line.replaceFirst('Error: ', '')}, enrichment skipped. "
                }
        }
        [pass, dir]
    }

    emit:
    result         = result_ch                                        // unchanged — published to 06.
    passed         = checked.filter { pass, _dir -> pass }.map { _pass, dir -> dir }   // for enrich (pass only)
    checkde_result = CHECKDE.out.checked.map { flag, _dir -> flag }   // QC receipt — published to 06.
    de_plot        = DE_PLOT.out.plot
    pca_plot       = PCA_PLOT.out.pca
}
