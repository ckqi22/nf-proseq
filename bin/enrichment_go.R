#!/usr/bin/env Rscript
# =============================================================================
# enrichment_go.R — GO enrichment analysis using clusterProfiler
#
# Reads edgeR DE results, runs enrichGO for each comparison (up/down separately),
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

argv <- arg_parser("GO enrichment analysis for PRO-seq DEGs")
argv <- add_argument(argv, "--de_results", help = "Diff_genes.xlsx from edger_de.R")
argv <- add_argument(argv, "--org_db",     help = "OrgDb package name", default = "org.Hs.eg.db")
argv <- add_argument(argv, "--output_dir", help = "Output directory for GO results")
argv <- parse_args(argv)

dir.create(argv$output_dir, recursive = TRUE, showWarnings = FALSE)

message("[enrichment_go] Reading DE results: ", argv$de_results)
de_data <- read.xlsx(argv$de_results, sheet = 1)
message("  ", nrow(de_data), " DEGs total")

if (nrow(de_data) == 0 || all(is.na(de_data$gene_id))) {
    message("[enrichment_go] No DEGs found; skipping GO enrichment.")
    quit(save = "no", status = 0)
}

# Get unique comparisons
comparisons <- unique(de_data$Comparison)
message("[enrichment_go] Comparisons: ", paste(comparisons, collapse = ", "))

# Load OrgDb
org_db <- get(argv$org_db, envir = asNamespace(argv$org_db))

for (comp in comparisons) {
    message("[enrichment_go] Processing: ", comp)

    comp_data <- de_data[de_data$Comparison == comp, ]

    # Up-regulated genes
    up_genes <- comp_data$gene_id[comp_data$Significant == "Up"]
    # Down-regulated genes
    down_genes <- comp_data$gene_id[comp_data$Significant == "Down"]

    for (direction in c("Up", "Down")) {
        genes <- if (direction == "Up") up_genes else down_genes
        if (length(genes) < 5) {
            message("  ", direction, ": too few genes (", length(genes), "), skipping")
            next
        }

        message("  ", direction, ": ", length(genes), " genes")

        # Convert gene symbols to ENTREZ IDs
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

        # GO enrichment
        for (ont in c("BP", "MF", "CC")) {
            ego <- tryCatch(
                enrichGO(
                    gene          = gene_entrez$ENTREZID,
                    OrgDb         = org_db,
                    ont           = ont,
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.05,
                    qvalueCutoff  = 0.2,
                    readable      = TRUE
                ),
                error = function(e) {
                    message("  WARNING: enrichGO (", ont, ") failed: ", e$message)
                    return(NULL)
                }
            )

            if (is.null(ego) || nrow(ego) == 0) {
                message("  ", ont, ": 0 enriched terms")
                next
            }

            message("  ", ont, ": ", nrow(ego), " enriched terms")

            prefix <- paste0(comp, "_", direction, "_GO_", ont)

            # Save table
            write.xlsx(as.data.frame(ego),
                       file = file.path(argv$output_dir, paste0(prefix, ".xlsx")),
                       rowNames = FALSE)

            # Dot plot
            pdf(file.path(argv$output_dir, paste0(prefix, "_dotplot.pdf")), width = 10, height = 8)
            if (nrow(ego) >= 1) {
                print(dotplot(ego, showCategory = 20, title = paste(comp, direction, ont)))
            }
            dev.off()

            # Bar plot
            pdf(file.path(argv$output_dir, paste0(prefix, "_barplot.pdf")), width = 10, height = 8)
            if (nrow(ego) >= 1) {
                print(barplot(ego, showCategory = 20, title = paste(comp, direction, ont)))
            }
            dev.off()
        }
    }
}

message("[enrichment_go] Done.")
