#!/usr/bin/env Rscript
# =============================================================================
# edger_de.R — Differential expression analysis using edgeR
#
# Reads gene body count matrix, runs edgeR exact test for each comparison,
# outputs Diff_genes.xlsx (significant DEGs) and All_Comparisons_genes.xlsx
# =============================================================================
suppressWarnings(suppressMessages({
    library(argparser)
    library(edgeR)
    library(openxlsx)
}))

argv <- arg_parser("edgeR differential expression analysis for PRO-seq data")
argv <- add_argument(argv, "--counts",       help = "Gene body count matrix (tab-separated, genes=rows, samples=columns)")
argv <- add_argument(argv, "--groups",       help = "YAML file mapping group names to sample lists")
argv <- add_argument(argv, "--comparisons",  help = "CSV file with columns: group1, group2")
argv <- add_argument(argv, "--fc_cutoff",    help = "Fold change cutoff",          default = 1.5)
argv <- add_argument(argv, "--pval_cutoff",  help = "P-value cutoff",              default = 0.05)
argv <- add_argument(argv, "--output_dir",   help = "Output directory")
argv <- parse_args(argv)

dir.create(argv$output_dir, recursive = TRUE, showWarnings = FALSE)

fc_cutoff   <- as.numeric(argv$fc_cutoff)
pval_cutoff <- as.numeric(argv$pval_cutoff)

# --- Read count matrix ---
message("[edger_de] Reading count matrix: ", argv$counts)
counts <- read.delim(argv$counts, header = TRUE, row.names = 1, stringsAsFactors = FALSE)
message("  ", nrow(counts), " genes x ", ncol(counts), " samples")

# --- Read groups YAML ---
message("[edger_de] Reading groups: ", argv$groups)
groups_raw <- read.delim(argv$groups, header = FALSE, stringsAsFactors = FALSE,
                         comment.char = "#", sep = ":", strip.white = TRUE)
group_names  <- trimws(groups_raw[, 1])
group_samples <- lapply(strsplit(trimws(groups_raw[, 2]), ","), trimws)
names(group_samples) <- group_names

# Build sample-to-group mapping
sample_groups <- character()
for (gp in names(group_samples)) {
    for (samp in group_samples[[gp]]) {
        sample_groups[samp] <- gp
    }
}
message("  Groups: ", paste(names(group_samples), "=", sapply(group_samples, length), "samples", collapse = ", "))

# --- Read comparisons ---
message("[edger_de] Reading comparisons: ", argv$comparisons)
comps <- read.delim(argv$comparisons, header = TRUE, stringsAsFactors = FALSE, sep = ",")
message("  ", nrow(comps), " comparison(s)")

# --- Match count columns to sample groups ---
count_samples <- colnames(counts)
# Handle featureCounts column naming: strip path, keep sample name
count_samples_clean <- gsub("\\.bam$", "", count_samples)
count_samples_clean <- basename(count_samples_clean)

group_vec <- sample_groups[count_samples_clean]
names(group_vec) <- count_samples

valid_samples <- !is.na(group_vec)
if (!all(valid_samples)) {
    message("[edger_de] WARNING: ", sum(!valid_samples), " samples not matched to groups: ",
            paste(count_samples[!valid_samples], collapse = ", "))
    group_vec <- group_vec[valid_samples]
    counts <- counts[, valid_samples, drop = FALSE]
}

group_vec <- factor(group_vec)
message("[edger_de] Group factor: ", paste(levels(group_vec), collapse = ", "))

# --- edgeR analysis ---
# Build DGEList
dge <- DGEList(counts = counts, group = group_vec)

# Filter low-expressed genes
keep <- filterByExpr(dge, group = group_vec)
dge <- dge[keep, , keep.lib.sizes = FALSE]
message("[edger_de] Filtered: ", sum(keep), " genes retained (", sum(!keep), " removed)")

# TMM normalization
dge <- calcNormFactors(dge, method = "TMM")

# Estimate dispersion
design <- model.matrix(~ 0 + group_vec)
colnames(design) <- levels(group_vec)
dge <- estimateDisp(dge, design)

# --- Run comparisons ---
all_results <- list()
diff_results <- list()

for (i in seq_len(nrow(comps))) {
    g1 <- comps$group1[i]
    g2 <- comps$group2[i]
    comp_name <- paste0(g1, "_vs_", g2)

    message("[edger_de] Comparing: ", comp_name)

    if (!g1 %in% levels(group_vec) || !g2 %in% levels(group_vec)) {
        message("  Skipping: group not found in data")
        next
    }

    # Exact test
    et <- exactTest(dge, pair = c(g2, g1))

    res <- topTags(et, n = Inf)$table
    res$gene_id <- rownames(res)
    res$Comparison <- comp_name

    # Mark significance
    res$Significant <- "NS"
    res$Significant[res$FDR < pval_cutoff & res$logFC > log2(fc_cutoff)]  <- "Up"
    res$Significant[res$FDR < pval_cutoff & res$logFC < -log2(fc_cutoff)] <- "Down"

    n_up   <- sum(res$Significant == "Up")
    n_down <- sum(res$Significant == "Down")
    message("  ", n_up, " up, ", n_down, " down (FC>", fc_cutoff, ", FDR<", pval_cutoff, ")")

    all_results[[comp_name]] <- res

    # Significant only
    sig <- res[res$Significant != "NS", ]
    if (nrow(sig) > 0) {
        diff_results[[comp_name]] <- sig
    }
}

# --- Write outputs ---
# All comparisons combined
all_combined <- do.call(rbind, all_results)
all_path <- file.path(argv$output_dir, "All_Comparisons_genes.xlsx")
write.xlsx(all_combined, file = all_path, sheetName = "All_Genes", rowNames = FALSE)
message("[edger_de] All genes: ", all_path)

# Significant genes only
diff_combined <- if (length(diff_results) > 0) do.call(rbind, diff_results) else data.frame()
diff_path <- file.path(argv$output_dir, "Diff_genes.xlsx")
if (nrow(diff_combined) > 0) {
    write.xlsx(diff_combined, file = diff_path, sheetName = "DEGs", rowNames = FALSE)
    message("[edger_de] Diff genes: ", diff_path, " (", nrow(diff_combined), " total)")
} else {
    # Write empty file with header
    empty <- data.frame(gene_id = character(), logFC = numeric(), logCPM = numeric(),
                        PValue = numeric(), FDR = numeric(), Comparison = character(),
                        Significant = character())
    write.xlsx(empty, file = diff_path, sheetName = "DEGs", rowNames = FALSE)
    message("[edger_de] No DEGs found; wrote empty Diff_genes.xlsx")
}

message("[edger_de] Done.")
