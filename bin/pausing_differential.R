#!/usr/bin/env Rscript
# =============================================================================
# pausing_differential.R — Differential pausing between groups (Wilcoxon)
#
# Input: PI table (from pausing_index.R) + groups YAML + comparisons CSV.
# Output: tab-separated differential table.
# =============================================================================
suppressWarnings(suppressMessages({
    library(argparser)
    library(yaml)
}))

# ── Parse args ──
argv <- arg_parser("Differential pausing analysis (Wilcoxon rank-sum)")
argv <- add_argument(argv, "--pi",          help = "PI table (from pausing_index.R)")
argv <- add_argument(argv, "--groups",      help = "Groups YAML (group_name -> [samples])")
argv <- add_argument(argv, "--comparisons", help = "Comparisons CSV (columns: group1, group2)")
argv <- add_argument(argv, "--output",      help = "Output differential table (tsv)")
argv <- parse_args(argv)

# ── Run comparisons ──
run_comparisons <- function(pi, groups, comps) {
    res <- data.frame(Comparison = character(), Group1 = character(), Group2 = character(),
                      Median_PI_G1 = numeric(), Median_PI_G2 = numeric(),
                      P_Value = numeric(), N_genes = integer(), stringsAsFactors = FALSE)

    for (i in seq_len(nrow(comps))) {
        g1 <- comps$group1[i]
        g2 <- comps$group2[i]
        s1 <- unlist(groups[[g1]])
        s2 <- unlist(groups[[g2]])
        if (is.null(s1) || is.null(s2)) next
        c1 <- paste0(s1, "_PI")
        c2 <- paste0(s2, "_PI")
        if (any(!c1 %in% colnames(pi)) || any(!c2 %in% colnames(pi))) next

        m1 <- rowMeans(pi[, c1, drop = FALSE], na.rm = TRUE)
        m2 <- rowMeans(pi[, c2, drop = FALSE], na.rm = TRUE)
        valid <- !is.na(m1) & !is.na(m2)
        if (sum(valid) < 2) next

        wt <- wilcox.test(m1[valid], m2[valid], paired = FALSE)
        res <- rbind(res, data.frame(
            Comparison   = paste0(g1, "_vs_", g2),
            Group1       = g1,
            Group2       = g2,
            Median_PI_G1 = median(m1[valid], na.rm = TRUE),
            Median_PI_G2 = median(m2[valid], na.rm = TRUE),
            P_Value      = wt$p.value,
            N_genes      = sum(valid),
            stringsAsFactors = FALSE))
    }
    res
}

# ── Main ──
main <- function(argv) {
    if (!file.exists(argv$pi)) stop("PI table not found: ", argv$pi)

    pi <- read.delim(argv$pi, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)

    groups <- NULL
    if (file.exists(argv$groups) && file.size(argv$groups) > 0) {
        groups <- yaml::read_yaml(argv$groups)
    }
    comps <- read.delim(argv$comparisons, header = TRUE, sep = ",", stringsAsFactors = FALSE)

    res <- if (is.null(groups) || nrow(comps) == 0) data.frame() else run_comparisons(pi, groups, comps)

    if (nrow(res) == 0) {
        res <- data.frame(Note = "No differential pausing: no comparisons configured",
                          stringsAsFactors = FALSE)
    }
    write.table(res, file = argv$output, sep = "\t", quote = FALSE, row.names = FALSE)
    message("[pausing_differential] -> ", argv$output)
}

main(argv)
