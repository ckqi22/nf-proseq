#!/usr/bin/env nextflow
//
// SUBWORKFLOW: diff_analysis
// Extract gene body columns → combined matrix → DE + enrichment + plots.
//

include { extract_region  } from '../modules/extract_region.nf'
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
    counts_files       // channel: per-sample proseq_counts.txt files
    comparisons_config // channel: val(list)
    config_ch          // channel: val(config)

    main:
    all_files = counts_files.collect()

    // Extract gene body matrix
    gb_mat = extract_region(all_files, "_gene_body")

    // Groups YAML
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

    // Comparisons CSV
    comp_csv = comparisons_config.map { comps ->
        def f = file("${workDir}/de_comparisons.csv")
        def compList = comps instanceof List ? comps : []
        f.text = "group1,group2\n"
        compList.each { c -> f.text += "${c[0]},${c[1]}\n" }
        return f
    }

    gb_matrix = gb_mat.matrix

    edger_de(gb_matrix, groups_yml, comp_csv, config_ch)
    enrichment_go(edger_de.out.diff_genes, config_ch)
    enrichment_kegg(edger_de.out.diff_genes, config_ch)
    gsea(gb_matrix, groups_yml, comp_csv, config_ch)
    pca_plot(gb_matrix, config_ch)
    heatmap(edger_de.out.diff_genes, gb_matrix, config_ch)
    scatter(gb_matrix, comp_csv, config_ch)
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
