#!/usr/bin/env nextflow
//
// SUBWORKFLOW: pause_analysis
// 暂停指数（计数长度归一化后的比值）：promoter / genebody。
//   输入为 pol2_count 产出的合并计数矩阵（promoter_matrix / genebody_matrix）。
//   分组来自 samplesheet。
//   1. 算 PI 表（pausing_index.R）
//   2. 按分组画 boxplot（pausing_boxplot.R）
//

include { pausing_index   } from '../modules/pause_analysis/pausing_index.nf'
include { pausing_boxplot } from '../modules/pause_analysis/pausing_boxplot.nf'

workflow pause_analysis {
    take:
    promoter_matrix   // channel: path(pol2_promoter.matrix.txt)  — gene_id, length, <samples>
    genebody_matrix   // channel: path(pol2_genebody.matrix.txt)
    groups_config     // channel: val(map) — group_name -> [samples] (from samplesheet)

    main:
    pausing_index(promoter_matrix, genebody_matrix)

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
