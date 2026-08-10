#!/usr/bin/env Rscript
#
# pausing_index.R - Compute PRO-seq pausing index for each gene
#
# Computes for each gene i:
#   PI_i = (TSS -50 to +300bp raw reads) / (GeneBody TSS+301bp to TES raw reads)
#
# No length normalization is applied — raw counts are used as-is per the PRO-seq plan.
# Genes with a gene body shorter than --min_gene_length (default 800 bp) are excluded.
#
# Usage:
#   Rscript pausing_index.R --tss_counts <tss.txt> --gb_counts <gb.txt> \
#       --groups <groups.yml> --comparisons <comparisons.csv> \
#       --output_dir <dir> [--spike_factors <factors.txt>] [--min_gene_length 800]
#

suppressPackageStartupMessages({
  library(argparse)
  library(data.table)
  library(openxlsx)
  library(ggplot2)
  library(yaml)
})

# ---------------------------------------------------------------------------
# Logging helpers
# ---------------------------------------------------------------------------
log_info <- function(msg) {
  cat(sprintf("[%s] INFO  %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg))
}

log_error <- function(msg) {
  cat(sprintf("[%s] ERROR %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg), file = stderr())
}

log_warn <- function(msg) {
  cat(sprintf("[%s] WARN  %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg))
}

# ---------------------------------------------------------------------------
# Read count matrix (tab-separated, first column = gene_id, then sample columns)
# ---------------------------------------------------------------------------
read_count_matrix <- function(path, label) {
  if (!file.exists(path)) {
    log_error(sprintf("File not found: %s (%s)", path, label))
    quit(status = 1)
  }
  log_info(sprintf("Reading %s: %s (%d rows x %d cols sampled)", label, path,
    length(readLines(path)) - 1, length(strsplit(readLines(path, n = 1), "\t")[[1]])))
  dt <- fread(path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  return(dt)
}

# ---------------------------------------------------------------------------
# Read spike-in factors
# ---------------------------------------------------------------------------
read_spike_factors <- function(path) {
  if (is.null(path) || !file.exists(path)) {
    log_warn("No spike factors file provided or not found; proceeding without spike normalization")
    return(NULL)
  }
  log_info(sprintf("Reading spike factors: %s", path))
  dt <- fread(path, header = TRUE, sep = "\t")
  log_info(sprintf("  %d samples with sigma/epsilon values", nrow(dt)))
  return(dt)
}

# ---------------------------------------------------------------------------
# Filter genes by minimum gene body length
# The gene_body file is expected to have columns: gene_id, length
# or we can infer length from the BED used to generate counts.
# Here we assume gb_counts has a 'length' column, or we read it from a supplemental file.
# ---------------------------------------------------------------------------
filter_by_gene_length <- function(gb_dt, min_length) {
  if ("length" %in% names(gb_dt)) {
    keep <- gb_dt$length >= min_length
    n_removed <- sum(!keep)
    if (n_removed > 0) {
      log_info(sprintf("Filtering genes: %d removed (gene body < %d bp), %d retained",
        n_removed, min_length, sum(keep)))
    }
    return(gb_dt[keep, ])
  } else {
    log_warn("No 'length' column found in gene body counts; skipping length filtering")
    return(gb_dt)
  }
}

# ---------------------------------------------------------------------------
# Compute pausing index: PI = TSS reads / GeneBody reads
# Matches samples by column name across the two matrices.
# ---------------------------------------------------------------------------
compute_pausing_index <- function(tss_dt, gb_dt, spike_factors, min_gene_length) {
  gene_col <- names(tss_dt)[1]
  if (names(gb_dt)[1] != gene_col) {
    log_error("TSS and GeneBody count matrices must have the same first column (gene ID)")
    quit(status = 1)
  }

  # Get gene IDs
  tss_genes <- tss_dt[[gene_col]]
  gb_genes  <- gb_dt[[gene_col]]
  common_genes <- intersect(tss_genes, gb_genes)
  log_info(sprintf("Common genes between TSS and GeneBody: %d", length(common_genes)))

  # Subset to common genes
  tss_dt <- tss_dt[tss_dt[[gene_col]] %in% common_genes, ]
  gb_dt  <- gb_dt[gb_dt[[gene_col]] %in% common_genes, ]

  # Identify shared sample columns
  tss_samples <- setdiff(names(tss_dt), gene_col)
  gb_samples  <- setdiff(names(gb_dt), gene_col)

  # If length column exists in gb, keep it separate
  length_col <- if ("length" %in% gb_samples) "length" else NULL
  gb_samples <- setdiff(gb_samples, "length")

  # Filter by min gene body length if length data exists
  if (!is.null(length_col)) {
    keep <- gb_dt[[length_col]] >= min_gene_length
    tss_dt <- tss_dt[keep, ]
    gb_dt  <- gb_dt[keep, ]
  }

  common_samples <- intersect(tss_samples, gb_samples)
  log_info(sprintf("Shared samples: %d", length(common_samples)))

  # Build results table
  genes <- tss_dt[[gene_col]]
  pi_dt <- data.table(gene_id = genes)

  for (samp in common_samples) {
    tss_vals <- as.numeric(tss_dt[[samp]])
    gb_vals  <- as.numeric(gb_dt[[samp]])

    promoter_reads_col <- paste0(samp, "_PromoterReads")
    genebody_reads_col <- paste0(samp, "_GeneBodyReads")
    pi_col             <- paste0(samp, "_PI")

    # Avoid division by zero
    pi_vals <- ifelse(gb_vals > 0, tss_vals / gb_vals, NA_real_)

    pi_dt[[promoter_reads_col]] <- tss_vals
    pi_dt[[genebody_reads_col]] <- gb_vals
    pi_dt[[pi_col]]             <- pi_vals

    log_info(sprintf("  %s: mean PI = %.4f (sd=%.4f)", samp,
      mean(pi_vals, na.rm = TRUE), sd(pi_vals, na.rm = TRUE)))
  }

  return(pi_dt)
}

# ---------------------------------------------------------------------------
# Read groups YAML. Expected format:
#   group_name:
#     - sample1
#     - sample2
# ---------------------------------------------------------------------------
read_groups <- function(path) {
  if (is.null(path) || !file.exists(path)) {
    log_warn("No groups file provided; skipping group comparisons")
    return(NULL)
  }
  log_info(sprintf("Reading groups: %s", path))
  groups <- read_yaml(path)
  return(groups)
}

# ---------------------------------------------------------------------------
# Read comparisons CSV. Expected columns: group1, group2
# ---------------------------------------------------------------------------
read_comparisons <- function(path) {
  if (is.null(path) || !file.exists(path)) {
    return(NULL)
  }
  log_info(sprintf("Reading comparisons: %s", path))
  comps <- fread(path, header = TRUE, sep = ",")
  return(comps)
}

# ---------------------------------------------------------------------------
# Run Wilcoxon rank-sum test for each comparison
# ---------------------------------------------------------------------------
run_comparisons <- function(pi_dt, groups, comparisons) {
  if (is.null(groups) || is.null(comparisons)) {
    log_info("Skipping differential pausing analysis (no groups/comparisons provided)")
    return(NULL)
  }

  diff_results <- data.table(
    Comparison  = character(),
    Group1      = character(),
    Group2      = character(),
    Median_PI_G1 = numeric(),
    Median_PI_G2 = numeric(),
    P_Value     = numeric(),
    N_genes     = integer()
  )

  for (i in seq_len(nrow(comparisons))) {
    g1_name <- comparisons$group1[i]
    g2_name <- comparisons$group2[i]

    g1_samples <- groups[[g1_name]]
    g2_samples <- groups[[g2_name]]

    if (is.null(g1_samples) || is.null(g2_samples)) {
      log_warn(sprintf("Skipping comparison '%s vs %s': group not found", g1_name, g2_name))
      next
    }

    # Pool PI values: for each sample, get the PI column
    g1_pi_cols <- paste0(g1_samples, "_PI")
    g2_pi_cols <- paste0(g2_samples, "_PI")

    missing_cols <- c(
      setdiff(g1_pi_cols, names(pi_dt)),
      setdiff(g2_pi_cols, names(pi_dt))
    )
    if (length(missing_cols) > 0) {
      log_warn(sprintf("Skipping comparison '%s vs %s': missing PI columns: %s",
        g1_name, g2_name, paste(missing_cols, collapse = ", ")))
      next
    }

    # Take mean PI per gene across replicates for each group
    g1_mean_pi <- rowMeans(pi_dt[, ..g1_pi_cols, with = FALSE], na.rm = TRUE)
    g2_mean_pi <- rowMeans(pi_dt[, ..g2_pi_cols, with = FALSE], na.rm = TRUE)

    # Remove NA pairs
    valid <- !is.na(g1_mean_pi) & !is.na(g2_mean_pi)
    g1_mean_pi <- g1_mean_pi[valid]
    g2_mean_pi <- g2_mean_pi[valid]

    wt <- wilcox.test(g1_mean_pi, g2_mean_pi, paired = FALSE, alternative = "two.sided")

    diff_results <- rbind(diff_results, data.table(
      Comparison   = paste0(g1_name, "_vs_", g2_name),
      Group1       = g1_name,
      Group2       = g2_name,
      Median_PI_G1 = median(g1_mean_pi, na.rm = TRUE),
      Median_PI_G2 = median(g2_mean_pi, na.rm = TRUE),
      P_Value      = wt$p.value,
      N_genes      = length(g1_mean_pi)
    ))

    log_info(sprintf("  %s vs %s: p = %.4e (Wilcoxon)", g1_name, g2_name, wt$p.value))
  }

  return(diff_results)
}

# ---------------------------------------------------------------------------
# Boxplot of pausing index by group
# ---------------------------------------------------------------------------
generate_boxplot <- function(pi_dt, groups, output_dir) {
  if (is.null(groups)) return(NULL)

  plot_data <- data.table(gene_id = character(), PI = numeric(), Group = character())

  for (grp_name in names(groups)) {
    samples <- groups[[grp_name]]
    pi_cols <- paste0(samples, "_PI")
    missing <- setdiff(pi_cols, names(pi_dt))
    if (length(missing) > 0) next

    # Mean PI per gene across group replicates
    grp_mean <- rowMeans(pi_dt[, ..pi_cols, with = FALSE], na.rm = TRUE)
    valid <- !is.na(grp_mean)

    plot_data <- rbind(plot_data, data.table(
      gene_id = pi_dt$gene_id[valid],
      PI      = grp_mean[valid],
      Group   = grp_name
    ))
  }

  if (nrow(plot_data) == 0) return(NULL)

  p <- ggplot(plot_data, aes(x = Group, y = log2(PI + 0.001), fill = Group)) +
    geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.3) +
    labs(
      title = "PRO-seq Pausing Index by Group",
      x     = "Group",
      y     = "log2(Pausing Index)"
    ) +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, face = "bold")
    ) +
    scale_fill_brewer(palette = "Set2")

  out_path <- file.path(output_dir, "Pausing_Index_Boxplot.pdf")
  ggsave(out_path, plot = p, width = 8, height = 6, dpi = 300)
  log_info(sprintf("Boxplot saved: %s", out_path))
  return(out_path)
}

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
main <- function() {
  parser <- ArgumentParser(
    description = "Compute PRO-seq pausing index for each gene"
  )
  parser$add_argument("--tss_counts",
    required = TRUE,
    help = "TSS window count matrix (-50 to +300 bp, tab-separated)")
  parser$add_argument("--gb_counts",
    required = TRUE,
    help = "Gene body count matrix (TSS+301 to TES, tab-separated)")
  parser$add_argument("--groups",
    default = NULL,
    help = "YAML file mapping group names to sample lists")
  parser$add_argument("--comparisons",
    default = NULL,
    help = "CSV file with columns: group1, group2")
  parser$add_argument("--output_dir",
    required = TRUE,
    help = "Directory for output files")
  parser$add_argument("--spike_factors",
    default = NULL,
    help = "Spike correction factors (tab-separated, from spike_correction.R)")
  parser$add_argument("--min_gene_length",
    type = "integer",
    default = 800L,
    help = "Minimum gene body length in bp (default: 800)")

  args <- parser$parse_args()

  # Create output directory
  dir.create(args$output_dir, showWarnings = FALSE, recursive = TRUE)

  log_info("=== PRO-seq Pausing Index Computation ===")

  # Read data
  tss_dt  <- read_count_matrix(args$tss_counts, "TSS counts")
  gb_dt   <- read_count_matrix(args$gb_counts,  "GeneBody counts")
  spike   <- read_spike_factors(args$spike_factors)
  groups  <- read_groups(args$groups)
  comps   <- read_comparisons(args$comparisons)

  # Compute pausing index
  pi_dt <- compute_pausing_index(tss_dt, gb_dt, spike, args$min_gene_length)
  log_info(sprintf("Pausing index computed for %d genes", nrow(pi_dt)))

  # 1. Write all genes output
  all_path <- file.path(args$output_dir, "Pausing_Index_All.xlsx")
  write.xlsx(pi_dt, file = all_path, sheetName = "Pausing_Index", rowNames = FALSE)
  log_info(sprintf("All genes output: %s", all_path))

  # 2. Boxplot
  generate_boxplot(pi_dt, groups, args$output_dir)

  # 3. Differential pausing
  diff_results <- run_comparisons(pi_dt, groups, comps)
  if (!is.null(diff_results) && nrow(diff_results) > 0) {
    diff_path <- file.path(args$output_dir, "Pausing_Differential.xlsx")
    write.xlsx(diff_results, file = diff_path, sheetName = "Differential_Pausing", rowNames = FALSE)
    log_info(sprintf("Differential pausing output: %s", diff_path))
  }

  log_info("Done.")
}

# ---------------------------------------------------------------------------
if (sys.nframe() == 0) {
  main()
}
