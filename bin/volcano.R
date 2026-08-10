#!/usr/bin/env Rscript
# =============================================================================
# volcano.R — Volcano plot for PRO-seq differential expression results
#
# Creates a volcano plot: -log10(p-value) vs log2(fold change).
# Significant genes colored: Up = red, Down = blue, NS = grey.
# =============================================================================
suppressWarnings(suppressMessages({
    library(argparser)
    library(ggplot2)
    library(openxlsx)
    library(ggrepel)
}))

argv <- arg_parser("Volcano plot for PRO-seq differential expression")
argv <- add_argument(argv, "--de_results",   help = "All_Comparisons_genes.xlsx from edger_de.R")
argv <- add_argument(argv, "--fc_cutoff",    help = "Fold change cutoff",       default = 1.5)
argv <- add_argument(argv, "--pval_cutoff",  help = "P-value cutoff (FDR)",     default = 0.05)
argv <- add_argument(argv, "--top_n_labels", help = "Number of top genes to label", default = 10)
argv <- add_argument(argv, "--output",       help = "Output PDF path")
argv <- parse_args(argv)

fc_cutoff   <- as.numeric(argv$fc_cutoff)
pval_cutoff <- as.numeric(argv$pval_cutoff)
top_n       <- as.integer(argv$top_n_labels)

# --- Read DE results ---
message("[volcano] Reading DE results: ", argv$de_results)
de_data <- read.xlsx(argv$de_results, sheet = 1)
message("  ", nrow(de_data), " genes")

if (nrow(de_data) == 0) {
    pdf(argv$output, width = 8, height = 6)
    plot(0, 0, type = "n", axes = FALSE, xlab = "", ylab = "")
    text(0, 0, "No data for volcano plot", cex = 1.5)
    dev.off()
    quit(save = "no", status = 0)
}

# Check required columns
required_cols <- c("gene_id", "logFC", "FDR")
for (col in required_cols) {
    if (!col %in% names(de_data)) {
        # Try PValue instead of FDR
        if (col == "FDR" && "PValue" %in% names(de_data)) {
            de_data$FDR <- de_data$PValue
            next
        }
        stop("Required column '", col, "' not found in DE results")
    }
}

# --- Classify genes ---
de_data$Significant <- "NS"
de_data$Significant[de_data$FDR < pval_cutoff & de_data$logFC > log2(fc_cutoff)]  <- "Up"
de_data$Significant[de_data$FDR < pval_cutoff & de_data$logFC < -log2(fc_cutoff)] <- "Down"

de_data$negLog10FDR <- -log10(de_data$FDR)
de_data$negLog10FDR[is.infinite(de_data$negLog10FDR)] <- max(
    de_data$negLog10FDR[is.finite(de_data$negLog10FDR)], na.rm = TRUE
)

n_up   <- sum(de_data$Significant == "Up")
n_down <- sum(de_data$Significant == "Down")
n_ns   <- sum(de_data$Significant == "NS")
message("[volcano] ", n_up, " up, ", n_down, " down, ", n_ns, " NS")

# --- Label top genes ---
sig_genes <- de_data[de_data$Significant != "NS", ]
sig_genes <- sig_genes[order(sig_genes$FDR), ]
label_genes <- head(sig_genes$gene_id, top_n)

de_data$Label <- ifelse(de_data$gene_id %in% label_genes, de_data$gene_id, "")

# --- Plot ---
color_map <- c("Up" = "#E64B35", "Down" = "#4DBBD5", "NS" = "#999999")

p <- ggplot(de_data, aes(x = logFC, y = negLog10FDR, color = Significant)) +
    geom_point(size = 0.8, alpha = 0.6) +
    geom_text_repel(aes(label = Label), size = 3, max.overlaps = 20,
                    box.padding = 0.5, point.padding = 0.3) +
    scale_color_manual(values = color_map) +
    geom_hline(yintercept = -log10(pval_cutoff), linetype = "dashed", color = "grey50") +
    geom_vline(xintercept =  log2(fc_cutoff),  linetype = "dashed", color = "grey50") +
    geom_vline(xintercept = -log2(fc_cutoff),  linetype = "dashed", color = "grey50") +
    labs(
        title    = "PRO-seq Differential Expression",
        subtitle = paste0("FC > ", fc_cutoff, ", FDR < ", pval_cutoff,
                          " | Up: ", n_up, " | Down: ", n_down),
        x        = expression(log[2]~Fold~Change),
        y        = expression(-log[10]~FDR),
        color    = ""
    ) +
    theme_bw(base_size = 12) +
    theme(
        plot.title    = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "bottom"
    )

pdf(argv$output, width = 8, height = 8)
print(p)
dev.off()

message("[volcano] Done: ", argv$output)
