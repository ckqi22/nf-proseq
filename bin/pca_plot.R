#!/usr/bin/env Rscript
# =============================================================================
# pca_plot.R — PCA plot from PRO-seq count matrix
#
# Performs PCA on log-CPM normalized counts, plots PC1 vs PC2 with sample labels.
# =============================================================================
suppressWarnings(suppressMessages({
    library(argparser)
    library(ggplot2)
    library(edgeR)
}))

argv <- arg_parser("PCA plot for PRO-seq data")
argv <- add_argument(argv, "--counts", help = "Count matrix (tab-separated, genes=rows, samples=columns)")
argv <- add_argument(argv, "--groups", help = "Optional YAML-like file mapping group names to sample lists")
argv <- add_argument(argv, "--output", help = "Output PDF path")
argv <- parse_args(argv)

# --- Read count matrix ---
message("[pca_plot] Reading count matrix: ", argv$counts)
counts <- read.delim(argv$counts, header = TRUE, row.names = 1, stringsAsFactors = FALSE)
message("  ", nrow(counts), " genes x ", ncol(counts), " samples")

if (ncol(counts) < 2) {
    message("[pca_plot] Only 1 sample; PCA requires >= 2 samples.")
    pdf(argv$output, width = 8, height = 6)
    plot(0, 0, type = "n", axes = FALSE, xlab = "", ylab = "")
    text(0, 0, "PCA requires >= 2 samples", cex = 1.5)
    dev.off()
    quit(save = "no", status = 0)
}

# --- Read groups if provided ---
group_labels <- NULL
if (!is.na(argv$groups) && argv$groups != "" && file.exists(argv$groups)) {
    message("[pca_plot] Reading groups: ", argv$groups)
    groups_raw <- read.delim(argv$groups, header = FALSE, stringsAsFactors = FALSE,
                             comment.char = "#", sep = ":", strip.white = TRUE)
    group_names  <- trimws(groups_raw[, 1])
    group_samples <- lapply(strsplit(trimws(groups_raw[, 2]), ","), trimws)
    names(group_samples) <- group_names

    sample_to_group <- character()
    for (gp in names(group_samples)) {
        for (samp in group_samples[[gp]]) {
            sample_to_group[samp] <- gp
        }
    }
    group_labels <- sample_to_group[colnames(counts)]
    message("  Groups: ", paste(unique(group_labels), collapse = ", "))
}

# --- Normalize and PCA ---
dge <- DGEList(counts = counts)
dge <- calcNormFactors(dge, method = "TMM")
log_cpm <- cpm(dge, log = TRUE, prior.count = 1)

# Remove rows with zero variance
row_vars <- apply(log_cpm, 1, var, na.rm = TRUE)
log_cpm <- log_cpm[row_vars > 0, , drop = FALSE]
message("[pca_plot] ", nrow(log_cpm), " genes with variance > 0")

pca <- prcomp(t(log_cpm), center = TRUE, scale. = TRUE)
var_explained <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)

plot_data <- data.frame(
    Sample = colnames(counts),
    PC1    = pca$x[, 1],
    PC2    = pca$x[, 2],
    Group  = if (!is.null(group_labels)) group_labels else colnames(counts),
    stringsAsFactors = FALSE
)

# --- Plot ---
p <- ggplot(plot_data, aes(x = PC1, y = PC2, color = Group, label = Sample)) +
    geom_point(size = 3) +
    geom_text_repel(size = 3.5, box.padding = 0.5) +
    labs(
        title = "PRO-seq PCA Plot",
        x     = paste0("PC1 (", var_explained[1], "%)"),
        y     = paste0("PC2 (", var_explained[2], "%)")
    ) +
    theme_bw(base_size = 12) +
    theme(
        plot.title      = element_text(hjust = 0.5, face = "bold"),
        legend.position = "bottom"
    )

if (is.null(group_labels)) {
    p <- p + theme(legend.position = "none")
}

pdf(argv$output, width = 8, height = 7)
print(p)
dev.off()

message("[pca_plot] Done: ", argv$output)
