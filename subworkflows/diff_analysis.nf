#!/usr/bin/env nextflow
//
// SUBWORKFLOW: diff_analysis
// Chains: merge counts -> edgeR DE -> GO/KEGG enrichment -> GSEA -> plots
//

include { edger_de        } from '../modules/edger_de.nf'
include { enrichment_go   } from '../modules/enrichment_go.nf'
include { enrichment_kegg } from '../modules/enrichment_kegg.nf'
include { gsea            } from '../modules/gsea.nf'
include { heatmap         } from '../modules/heatmap.nf'
include { scatter         } from '../modules/scatter.nf'
include { volcano         } from '../modules/volcano.nf'
include { pca_plot        } from '../modules/pca_plot.nf'

workflow diff_analysis {
    take:
    gene_body_counts   // channel: per-sample gene body count files
    comparisons_config // channel: val(list)
    config_ch          // channel: val(config)

    main:
    // Collect per-sample files
    gb_files = gene_body_counts.collect()

    // Merge into combined matrix
    merged_counts = gb_files.map { files ->
        def out = file("${workDir}/gb_merged_matrix.txt")
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

    // Generate groups YAML
    groups_yml = Channel.value(params.group ?: [:]).map { groups ->
        def lines = []
        groups.each { groupName, sampleList ->
            lines << "${groupName}: ${sampleList}"
        }
        if (lines.isEmpty()) { lines << "all: unknown" }
        def f = file("${workDir}/de_groups.yml")
        f.text = lines.join('\n')
        return f
    }

    // Generate comparisons CSV
    comp_csv = comparisons_config.map { comps ->
        def f = file("${workDir}/de_comparisons.csv")
        def compList = comps instanceof List ? comps : []
        f.text = "group1,group2\n"
        compList.each { c -> f.text += "${c[0]},${c[1]}\n" }
        return f
    }

    // edgeR DE
    edger_de(merged_counts, groups_yml, comp_csv, config_ch)

    // GO / KEGG enrichment
    enrichment_go(edger_de.out.diff_genes, config_ch)
    enrichment_kegg(edger_de.out.diff_genes, config_ch)

    // GSEA
    gsea(merged_counts, groups_yml, comp_csv, config_ch)

    // Plots
    pca_plot(merged_counts, config_ch)
    heatmap(edger_de.out.diff_genes, merged_counts, config_ch)
    scatter(merged_counts, comp_csv, config_ch)
    volcano(edger_de.out.all_comparisons, config_ch)

    emit:
    diff_genes      = edger_de.out.diff_genes
    all_comparisons = edger_de.out.all_comparisons
    go_results      = enrichment_go.out.go_results
    kegg_results    = enrichment_kegg.out.kegg_results
    gsea_results    = gsea.out.gsea_results
    heatmap_pdf     = heatmap.out.heatmap_pdf
    scatter_pdf     = scatter.out.scatter_pdf
    volcano_pdf     = volcano.out.volcano_pdf
    pca_pdf         = pca_plot.out.pca_pdf
}
