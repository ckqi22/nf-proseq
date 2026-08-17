#!/usr/bin/env nextflow
//
// SUBWORKFLOW: pause_analysis
// Pausing index (length-normalized): TSS / gene body.
//   Inputs are the merged count matrices produced by `quantification`.
//   Groups come from the samplesheet.
//   1. compute PI table (pausing_index.R)
//   2. boxplot by group (pausing_boxplot.R)
//
//   NOTE: differential analysis no longer lives here — it moved to
//   subworkflows/diff.nf (DESeq2 on gene body counts).
//

include { pausing_index   } from '../modules/pause_analysis/pausing_index.nf'
include { pausing_boxplot } from '../modules/pause_analysis/pausing_boxplot.nf'

workflow pause_analysis {
    take:
    tss_matrix        // channel: path(tss.matrix.txt)  — gene_id, length, <samples>
    genebody_matrix   // channel: path(genebody.matrix.txt)
    groups_config     // channel: val(map) — group_name -> [samples] (from samplesheet)

    main:
    pausing_index(tss_matrix, genebody_matrix)

    // Build the groups YAML as a plain string (a val channel).
    // NOTE: never call file() / write files inside a .map closure — that triggers
    // a DataflowExpression.invokeMethod StackOverflowError. The actual file is
    // written inside each process' script block instead.
    groups_yml = groups_config.map { groups ->
        def lines = []
        groups.each { groupName, samples ->
            lines << "${groupName}:"
            samples.each { s -> lines << "  - ${s}" }
        }
        return lines.join('\n')
    }

    pausing_boxplot(pausing_index.out.pi, groups_yml)

    emit:
    pi_all     = pausing_index.out.pi
    pi_boxplot = pausing_boxplot.out.pdf
}
