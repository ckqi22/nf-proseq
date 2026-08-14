#!/usr/bin/env Rscript
# =============================================================================
# pausing_boxplot.R — Boxplot of per-gene pausing index by group
#
# Input: PI table (from pausing_index.R) + groups YAML (group_name -> [samples]).
# Output: a single PDF.
# =============================================================================
suppressWarnings(suppressMessages({
    library(argparser)
    library(ggplot2)
    library(yaml)
}))

# ── Parse args ──
argv <- arg_parser("Boxplot of PRO-seq pausing index by group")
argv <- add_argument(argv, "--pi",     help = "PI table (from pausing_index.R)")
argv <- add_argument(argv, "--groups", help = "Groups YAML (group_name -> [samples])")
argv <- add_argument(argv, "--output", help = "Output PDF")
argv <- parse_args(argv)

# ── Build plot data ──
build_plot_data <- function(pi, groups) {
    plot_data <- data.frame(gene_id = character(), PI = numeric(), Group = character(), stringsAsFactors = FALSE)
    for (grp in names(groups)) {
        samples <- unlist(groups[[grp]])
        pi_cols <- paste0(samples, "_PI")
        if (any(!pi_cols %in% colnames(pi))) next
        grp_mean <- rowMeans(pi[, pi_cols, drop = FALSE], na.rm = TRUE)
        valid    <- !is.na(grp_mean)
        plot_data <- rbind(plot_data, data.frame(
            gene_id = pi$gene_id[valid], PI = grp_mean[valid], Group = grp,
            stringsAsFactors = FALSE))
    }
    plot_data
}

# ── Main ──
main <- function(argv) {
    if (!file.exists(argv$pi)) stop("PI table not found: ", argv$pi)

    pi <- read.delim(argv$pi, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)

    groups <- NULL
    if (file.exists(argv$groups) && file.size(argv$groups) > 0) {
        groups <- yaml::read_yaml(argv$groups)
    }

    plot_data <- if (is.null(groups)) data.frame() else build_plot_data(pi, groups)

    if (nrow(plot_data) == 0) {
        pdf(argv$output, width = 8, height = 6)
        plot.new()
        text(0.5, 0.5, "No groups configured")
        dev.off()
    } else {
        p <- ggplot(plot_data, aes(x = Group, y = log2(PI + 0.001), fill = Group)) +
            geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.3) +
            labs(title = "PRO-seq Pausing Index by Group",
                 x = "Group", y = "log2(Pausing Index)") +
            theme_bw(base_size = 12) +
            theme(legend.position = "none",
                  plot.title = element_text(hjust = 0.5, face = "bold")) +
            scale_fill_brewer(palette = "Set2")
        ggsave(argv$output, plot = p, width = 8, height = 6, dpi = 300)
    }
    message("[pausing_boxplot] -> ", argv$output)
}

main(argv)
