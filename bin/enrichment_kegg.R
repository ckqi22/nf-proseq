#!/usr/bin/env Rscript
# =============================================================================
# enrichment_kegg.R — KEGG pathway enrichment using clusterProfiler
#
# Reads edgeR DE results, runs enrichKEGG for each comparison (up/down separately),
# generates dot plots and bar plots.
# =============================================================================
suppressWarnings(suppressMessages({
    library(argparser)
    library(clusterProfiler)
    library(org.Hs.eg.db)
    library(openxlsx)
    library(ggplot2)
    library(enrichplot)
}))

argv <- arg_parser("KEGG pathway enrichment analysis for PRO-seq DEGs")
argv <- add_argument(argv, "--de_results", help = "Diff_genes.xlsx from edger_de.R")
argv <- add_argument(argv, "--org_db",     help = "OrgDb package name", default = "org.Hs.eg.db")
argv <- add_argument(argv, "--organism",   help = "KEGG organism code", default = "hsa")
argv <- add_argument(argv, "--output_dir", help = "Output directory for KEGG results")
argv <- parse_args(argv)

dir.create(argv$output_dir, recursive = TRUE, showWarnings = FALSE)

message("[enrichment_kegg] Reading DE results: ", argv$de_results)
de_data <- read.xlsx(argv$de_results, sheet = 1)
message("  ", nrow(de_data), " DEGs total")

if (nrow(de_data) == 0 || all(is.na(de_data$gene_id))) {
    message("[enrichment_kegg] No DEGs found; skipping KEGG enrichment.")
    quit(save = "no", status = 0)
}

comparisons <- unique(de_data$Comparison)
message("[enrichment_kegg] Comparisons: ", paste(comparisons, collapse = ", "))

org_db <- get(argv$org_db, envir = asNamespace(argv$org_db))

for (comp in comparisons) {
    message("[enrichment_kegg] Processing: ", comp)

    comp_data <- de_data[de_data$Comparison == comp, ]

    for (direction in c("Up", "Down")) {
        genes <- if (direction == "Up") {
            comp_data$gene_id[comp_data$Significant == "Up"]
        } else {
            comp_data$gene_id[comp_data$Significant == "Down"]
        }

        if (length(genes) < 5) {
            message("  ", direction, ": too few genes (", length(genes), "), skipping")
            next
        }

        message("  ", direction, ": ", length(genes), " genes")

        # Convert symbols to ENTREZ IDs
        gene_entrez <- tryCatch(
            bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org_db),
            error = function(e) {
                message("  WARNING: gene ID conversion failed: ", e$message)
                return(NULL)
            }
        )

        if (is.null(gene_entrez) || nrow(gene_entrez) == 0) {
            message("  No valid ENTREZ IDs; skipping")
            next
        }

        # KEGG enrichment
        ekegg <- tryCatch(
            enrichKEGG(
                gene          = gene_entrez$ENTREZID,
                organism      = argv$organism,
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.2
            ),
            error = function(e) {
                message("  WARNING: enrichKEGG failed: ", e$message)
                return(NULL)
            }
        )

        if (is.null(ekegg) || nrow(ekegg) == 0) {
            message("  KEGG: 0 enriched pathways")
            next
        }

        message("  KEGG: ", nrow(ekegg), " enriched pathways")

        prefix <- paste0(comp, "_", direction, "_KEGG")

        # Save table
        write.xlsx(as.data.frame(ekegg),
                   file = file.path(argv$output_dir, paste0(prefix, ".xlsx")),
                   rowNames = FALSE)

        # Dot plot
        pdf(file.path(argv$output_dir, paste0(prefix, "_dotplot.pdf")), width = 10, height = 8)
        if (nrow(ekegg) >= 1) {
            print(dotplot(ekegg, showCategory = 20, title = paste(comp, direction, "KEGG")))
        }
        dev.off()

        # Bar plot
        pdf(file.path(argv$output_dir, paste0(prefix, "_barplot.pdf")), width = 10, height = 8)
        if (nrow(ekegg) >= 1) {
            print(barplot(ekegg, showCategory = 20, title = paste(comp, direction, "KEGG")))
        }
        dev.off()
    }
}

message("[enrichment_kegg] Done.")
