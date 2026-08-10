#!/usr/bin/env Rscript
# =============================================================================
# scatter.R — Sample correlation scatter plots for PRO-seq data
#
# Computes pairwise correlations between samples and generates scatter plots.
# =============================================================================
suppressWarnings(suppressMessages({
    library(argparser)
    library(ggplot2)
    library(GGally)
}))

argv <- arg_parser("Sample correlation scatter plots for PRO-seq data")
argv <- add_argument(argv, "--rpkm",   help = "Normalized expression matrix (RPKM/CPM, tab-separated)")
argv <- add_argument(argv, "--output", help = "Output PDF path")
argv <- parse_args(argv)

# --- Read expression matrix ---
message("[scatter] Reading expression matrix: ", argv$rpkm)
rpkm <- read.delim(argv$rpkm, header = TRUE, row.names = 1, stringsAsFactors = FALSE)
message("  ", nrow(rpkm), " genes x ", ncol(rpkm), " samples")

if (ncol(rpkm) < 2) {
    message("[scatter] Only 1 sample; cannot create scatter plot.")
    pdf(argv$output, width = 8, height = 6)
    plot(0, 0, type = "n", axes = FALSE, xlab = "", ylab = "")
    text(0, 0, "Only 1 sample — scatter plot requires >= 2 samples", cex = 1.2)
    dev.off()
    quit(save = "no", status = 0)
}

# --- Log2 transform ---
mat <- log2(as.matrix(rpkm) + 1)

# --- Create pairwise scatter matrix ---
message("[scatter] Generating scatter matrix for ", ncol(mat), " samples...")

pdf(argv$output, width = max(10, ncol(mat) * 3), height = max(8, ncol(mat) * 3))

if (ncol(mat) <= 6) {
    # Use GGally for nice scatter matrix
    plot_data <- as.data.frame(mat)
    p <- ggpairs(plot_data,
        lower = list(continuous = wrap("points", size = 0.3, alpha = 0.3)),
        upper = list(continuous = wrap("cor", size = 4)),
        diag  = list(continuous = wrap("densityDiag", alpha = 0.5)),
        title = "PRO-seq Sample Correlation"
    )
    print(p)
} else {
    # Too many samples — create a correlation heatmap
    cor_mat <- cor(mat, method = "pearson", use = "pairwise.complete.obs")
    cor_df <- reshape2::melt(cor_mat)
    names(cor_df) <- c("Sample1", "Sample2", "Correlation")

    p <- ggplot(cor_df, aes(x = Sample1, y = Sample2, fill = Correlation)) +
        geom_tile() +
        geom_text(aes(label = round(Correlation, 2)), size = 3) +
        scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0.5,
                             limits = c(0, 1)) +
        labs(title = "PRO-seq Sample Correlation Matrix",
             x = "", y = "") +
        theme_minimal(base_size = 11) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              panel.grid = element_blank())
    print(p)
}

dev.off()
message("[scatter] Done: ", argv$output)
