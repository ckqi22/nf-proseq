#!/usr/bin/env nextflow
//
// SUBWORKFLOW: pause_analysis
// From per-sample merged count files (3 regions), extract TSS and gene body
// columns, build combined matrices, and compute pausing index.
//

include { pausing_index } from '../modules/pausing_index.nf'

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

    // Extract TSS and gene body matrices from merged files
    matrices = all_files.map { files ->
        def tss_out = file("${workDir}/tss_merged.txt")
        def gb_out  = file("${workDir}/gb_merged.txt")
        def rf = files.collect { "'${it}'" }.join(', ')
        """
        ${params.r} -e "
            files <- c(${rf})
            first <- read.delim(files[1], header=TRUE, stringsAsFactors=FALSE)

            # TSS matrix: gene_id + _tss columns
            tss_col  <- 'gene_id'
            tss_cols <- grep('_tss\\$', names(first), value=TRUE)
            tss_mat  <- first[, c(tss_col, tss_cols), drop=FALSE]

            # Gene body matrix: gene_id + _gene_body columns
            gb_col   <- 'gene_id'
            gb_cols  <- grep('_gene_body\\$', names(first), value=TRUE)
            gb_mat   <- first[, c(gb_col, gb_cols), drop=FALSE]

            # Fill from remaining files
            for (f in files[-1]) {
                dt <- read.delim(f, header=TRUE, stringsAsFactors=FALSE)
                tc <- grep('_tss\\$', names(dt), value=TRUE)
                gc <- grep('_gene_body\\$', names(dt), value=TRUE)
                tss_mat[[tc[1]]] <- dt[[tc[1]]][match(tss_mat\\$gene_id, dt\\$gene_id)]
                gb_mat[[gc[1]]]   <- dt[[gc[1]]][match(gb_mat\\$gene_id, dt\\$gene_id)]
            }
            tss_mat[is.na(tss_mat)] <- 0
            gb_mat[is.na(gb_mat)]   <- 0

            write.table(tss_mat, file='${tss_out}', sep='\\t', quote=FALSE, row.names=FALSE)
            write.table(gb_mat,  file='${gb_out}',  sep='\\t', quote=FALSE, row.names=FALSE)
            message('[pause_analysis] TSS: ', nrow(tss_mat), ' genes x ', ncol(tss_mat)-1, ' samples')
            message('[pause_analysis] GB:  ', nrow(gb_mat),  ' genes x ', ncol(gb_mat)-1,  ' samples')
        "
        """
        return [tss_out, gb_out]
    }

    merged_tss = matrices.map { m -> m[0] }
    merged_gb  = matrices.map { m -> m[1] }

    pausing_index(merged_tss, merged_gb, groups_yml, Channel.empty())

    emit:
    pi_all     = pausing_index.out.pi_all
    pi_boxplot = pausing_index.out.pi_boxplot
    pi_diff    = pausing_index.out.pi_diff
}
