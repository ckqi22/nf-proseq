#!/usr/bin/env nextflow
//
// SUBWORKFLOW: pause_analysis
// Extract TSS and gene body columns → combined matrices → pausing index.
//

include { extract_region as extract_tss_region} from '../modules/extract_region.nf'
include { extract_region as extract_gb_region} from '../modules/extract_region.nf'
include { pausing_index   } from '../modules/pausing_index.nf'

workflow pause_analysis {
    take:
    counts_files       // channel: per-sample proseq_counts.txt files
    groups_config      // channel: val(map)

    main:
    all_files = counts_files.collect()

    // Generate groups YAML
    groups_yml = groups_config.map { groups ->
        def lines = []
        groups.each { groupName, sampleList ->
            lines << "${groupName}: ${sampleList}"
        }
        if (lines.isEmpty()) { lines << "all: unknown" }
        def f = file("${workDir}/pause_groups.yml")
        f.text = lines.join('\n')
        return f
    }

    // Extract TSS and gene body matrices
    tss_mat = extract_tss_region(all_files, "_tss")
    gb_mat  = extract_gb_region(all_files, "_gene_body")

    pausing_index(tss_mat.out.matrix, gb_mat.out.matrix, groups_yml, Channel.empty())

    emit:
    pi_all     = pausing_index.out.pi_all
    pi_boxplot = pausing_index.out.pi_boxplot
    pi_diff    = pausing_index.out.pi_diff
}
