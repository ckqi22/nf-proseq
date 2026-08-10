#!/usr/bin/env Rscript
# =============================================================================
# heatmap.R — Heatmap of differentially expressed genes using ComplexHeatmap
#
# Subsets expression matrix to top DEGs and generates a clustered heatmap.
# =============================================================================
suppressWarnings(suppressMessages({
    library(argparser)
    library(ComplexHeatmap)
    library(circlize)
    library(openxlsx)
}))

argv <- arg_parser("Heatmap of differentially expressed genes for PRO-seq data")
argv <- add_argument(argv, "--de_results", help = "Diff_genes.xlsx from edger_de.R")
argv <- add_argument(argv, "--rpkm",       help = "Normalized expression matrix (RPKM/CPM, tab-separated)")
argv <- add_argument(argv, "--top_n",      help = "Maximum number of genes to show", default = 500)
argv <- add_argument(argv, "--output",     help = "Output PDF path")
argv <- parse_args(argv)

top_n <- as.integer(argv$top_n)

# --- Read DE results ---
message("[heatmap] Reading DE results: ", argv$de_results)
de_data <- read.xlsx(argv$de_results, sheet = 1)
message("  ", nrow(de_data), " DEGs")

if (nrow(de_data) == 0 || all(is.na(de_data$gene_id))) {
    message("[heatmap] No DEGs found; skipping heatmap.")
    pdf(argv$output, width = 8, height = 6)
    plot(0, 0, type = "n", axes = FALSE, xlab = "", ylab = "")
    text(0, 0, "No DEGs to display", cex = 1.5)
    dev.off()
    quit(save = "no", status = 0)
}

# --- Read expression matrix ---
message("[heatmap] Reading expression matrix: ", argv$rpkm)
rpkm <- read.delim(argv$rpkm, header = TRUE, row.names = 1, stringsAsFactors = FALSE)
message("  ", nrow(rpkm), " genes x ", ncol(rpkm), " samples")

# --- Get top DEGs by significance ---
de_sorted <- de_data[order(de_data$PValue), ]
top_genes <- unique(de_sorted$gene_id)[1:min(top_n, length(unique(de_sorted$gene_id)))]
top_genes <- top_genes[!is.na(top_genes)]

# Match to expression matrix
present_genes <- intersect(top_genes, rownames(rpkm))
message("[heatmap] ", length(present_genes), " of ", length(top_genes), " top DEGs found in expression matrix")

if (length(present_genes) < 3) {
    message("[heatmap] Too few genes for heatmap; writing placeholder.")
    pdf(argv$output, width = 8, height = 6)
    plot(0, 0, type = "n", axes = FALSE, xlab = "", ylab = "")
    text(0, 0, paste("Only", length(present_genes), "genes — too few for heatmap"), cex = 1.5)
    dev.off()
    quit(save = "no", status = 0)
}

# Subset and z-score normalize by row
mat <- as.matrix(rpkm[present_genes, , drop = FALSE])
mat_log <- log2(mat + 1)
mat_z <- t(scale(t(mat_log)))

# Cap extreme values
mat_z[mat_z > 3]  <- 3
mat_z[mat_z < -3] <- -3

# Color mapping
col_fun <- colorRamp2(
    breaks = c(-3, 0, 3),
    colors = c("blue", "white", "red")
)

# Draw heatmap
message("[heatmap] Generating heatmap...")
pdf(argv$output, width = max(10, ncol(mat) * 1.5), height = max(8, length(present_genes) * 0.04))

ht <- Heatmap(mat_z,
    name                = "Z-score",
    col                 = col_fun,
    cluster_rows        = TRUE,
    cluster_columns     = TRUE,
    show_row_names      = length(present_genes) <= 100,
    show_column_names   = TRUE,
    column_title        = "PRO-seq DEGs Heatmap",
    row_title           = paste0(length(present_genes), " genes"),
    use_raster           = TRUE,
    heatmap_legend_param = list(direction = "horizontal")
)
draw(ht, heatmap_legend_side = "bottom")

dev.off()
message("[heatmap] Done: ", argv$output)
