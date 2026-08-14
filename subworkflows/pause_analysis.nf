#!/usr/bin/env nextflow
//
// SUBWORKFLOW: pause_analysis
// Pausing index (length-normalized): TSS / gene body.
//   1. convert featureCounts -> clean matrices
//   2. compute PI table (pausing_index.R)
//   3. boxplot by group + differential pausing (separate scripts)
//

include { FEATURECOUNTS_TO_MATRIX as TSS_MATRIX } from '../modules/featurecounts_to_matrix.nf'
include { FEATURECOUNTS_TO_MATRIX as GB_MATRIX }  from '../modules/featurecounts_to_matrix.nf'
include { pausing_index        } from '../modules/pausing_index.nf'
include { pausing_boxplot      } from '../modules/pausing_boxplot.nf'
include { pausing_differential } from '../modules/pausing_differential.nf'

workflow pause_analysis {
    take:
    tss_counts        // channel: path(tss.featureCounts.txt)
    genebody_counts   // channel: path(genebody.featureCounts.txt)
    groups_config     // channel: val(map) — group_name -> "sample1, sample2"

    main:
    TSS_MATRIX(tss_counts)
    GB_MATRIX(genebody_counts)

    pausing_index(TSS_MATRIX.out.matrix, GB_MATRIX.out.matrix)

    // Build the groups YAML / comparisons CSV as plain strings (a val channel).
    // NOTE: never call file() / write files inside a .map closure — that triggers
    // a DataflowExpression.invokeMethod StackOverflowError. The actual file is
    // written inside each process' script block instead.
    groups_yml = groups_config.map { groups ->
        def lines = []
        groups.each { groupName, sampleList ->
            def samples = sampleList.toString().split(',').collect { it.trim() }.findAll { it }
            lines << "${groupName}:"
            samples.each { lines << "  - ${it}" }
        }
        return lines.join('\n')
    }

    // TODO(comparisons): build group1,group2 CSV from params.compared_groups.
    // Header-only CSV for now -> no differential comparisons.
    comparisons_csv = groups_config.map { _ -> "group1,group2\n" }

    // Multicast the single PI table and groups to the two downstream analyses.
    pausing_index.out.pi.into { pi_all_ch; pi_box_ch; pi_diff_ch }
    groups_yml.into { groups_box_ch; groups_diff_ch }

    pausing_boxplot(pi_box_ch, groups_box_ch)
    pausing_differential(pi_diff_ch, groups_diff_ch, comparisons_csv)

    emit:
    pi_all     = pi_all_ch
    pi_boxplot = pausing_boxplot.out.pdf
    pi_diff    = pausing_differential.out.diff
}
