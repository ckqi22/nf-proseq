#!/usr/bin/env nextflow
//
// SUBWORKFLOW: pause_analysis
// Merges per-sample TSS and gene body counts, then computes pausing index.
//

include { pausing_index } from '../modules/pausing_index.nf'

workflow pause_analysis {
    take:
    tss_counts          // channel: per-sample TSS count files
    gene_body_counts    // channel: per-sample gene body count files
    groups_config       // channel: val(map) — group_name -> [sample1, sample2, ...]

    main:
    // Collect per-sample files
    tss_files = tss_counts.collect()
    gb_files  = gene_body_counts.collect()

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

    // Merge per-sample count files into combined matrices (R inline)
    merged_tss = tss_files.map { files ->
        def out = file("${workDir}/tss_merged.txt")
        def rf = files.collect { "'${it}'" }.join(', ')
        """
        ${params.r} -e "
            files <- c(${rf})
            first <- read.delim(files[1], header=TRUE, stringsAsFactors=FALSE)
            gene_col <- names(first)[1]
            result <- first[, gene_col, drop=FALSE]
            if ('Length' %in% names(first)) result\\$Length <- first\\$Length
            for (f in files) {
                nm <- sub('\\\\.tss_counts\\\\.txt\\$', '', basename(f))
                dt <- read.delim(f, header=TRUE, stringsAsFactors=FALSE)
                count_cols <- grep('_Total\\$', names(dt), value=TRUE)
                if (length(count_cols) == 0) count_cols <- tail(names(dt), 1)
                result[[nm]] <- dt[[count_cols[1]]]
            }
            write.table(result, file='${out}', sep='\\t', quote=FALSE, row.names=FALSE)
        "
        """
        return out
    }

    merged_gb = gb_files.map { files ->
        def out = file("${workDir}/gb_merged.txt")
        def rf = files.collect { "'${it}'" }.join(', ')
        """
        ${params.r} -e "
            files <- c(${rf})
            first <- read.delim(files[1], header=TRUE, stringsAsFactors=FALSE)
            gene_col <- names(first)[1]
            result <- first[, gene_col, drop=FALSE]
            if ('Length' %in% names(first)) result\\$Length <- first\\$Length
            for (f in files) {
                nm <- sub('\\\\.gene_body_counts\\\\.txt\\$', '', basename(f))
                dt <- read.delim(f, header=TRUE, stringsAsFactors=FALSE)
                count_cols <- grep('_Total\\$', names(dt), value=TRUE)
                if (length(count_cols) == 0) count_cols <- tail(names(dt), 1)
                result[[nm]] <- dt[[count_cols[1]]]
            }
            write.table(result, file='${out}', sep='\\t', quote=FALSE, row.names=FALSE)
        "
        """
        return out
    }

    // Run pausing index (no spike-in)
    pausing_index(merged_tss, merged_gb, groups_yml, Channel.empty())

    emit:
    pi_all     = pausing_index.out.pi_all
    pi_boxplot = pausing_index.out.pi_boxplot
    pi_diff    = pausing_index.out.pi_diff
}
