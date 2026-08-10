#!/usr/bin/env Rscript
#
# spike_correction.R - Compute spike-in normalization factors for PRO-seq data
#
# For each sample, computes:
#   sigma_j (sequencing depth factor) = spike_reads / endogenous_reads
#   epsilon_j (cross-contamination rate) if Input/Nascent paired data exists
#
# Usage:
#   Rscript spike_correction.R --endogenous <counts.txt> --spike <spike_counts.txt> --output <factors.txt>
#   [--epsilon_input <input_counts.txt> --epsilon_nascent <nascent_counts.txt>]
#

suppressPackageStartupMessages({
  library(argparse)
  library(data.table)
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
# Read count matrix (tab-separated, genes=rows, samples=columns)
# The first column should be gene_id; remaining columns are sample counts.
# ---------------------------------------------------------------------------
read_counts <- function(path, label) {
  if (!file.exists(path)) {
    log_error(sprintf("File not found: %s (%s)", path, label))
    quit(status = 1)
  }
  log_info(sprintf("Reading %s counts from: %s", label, path))
  dt <- fread(path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  # Ensure the first column is the gene/feature ID
  gene_col <- names(dt)[1]
  log_info(sprintf("  %d genes x %d samples", nrow(dt), ncol(dt) - 1))
  return(list(data = dt, gene_col = gene_col))
}

# ---------------------------------------------------------------------------
# Compute sigma_j = spike_reads / endogenous_reads for each sample
# ---------------------------------------------------------------------------
compute_sigma <- function(endo_counts, spike_counts) {
  endo_dt   <- endo_counts$data
  spike_dt  <- spike_counts$data
  gene_col  <- endo_counts$gene_col

  # Identify sample columns (everything except the gene column)
  sample_cols <- setdiff(names(endo_dt), gene_col)

  sigma <- data.table(sample = character(), sigma = numeric())

  for (samp in sample_cols) {
    if (!(samp %in% names(spike_dt))) {
      log_warn(sprintf("Sample '%s' not found in spike counts; setting sigma = NA", samp))
      sigma <- rbind(sigma, data.table(sample = samp, sigma = NA_real_))
      next
    }

    endo_total  <- sum(endo_dt[[samp]], na.rm = TRUE)
    spike_total <- sum(spike_dt[[samp]], na.rm = TRUE)

    if (endo_total == 0) {
      log_warn(sprintf("Sample '%s' has zero endogenous reads; setting sigma = NA", samp))
      sigma <- rbind(sigma, data.table(sample = samp, sigma = NA_real_))
      next
    }

    s <- spike_total / endo_total
    log_info(sprintf("  %s: sigma = %.6f (spike=%d, endo=%d)", samp, s, spike_total, endo_total))
    sigma <- rbind(sigma, data.table(sample = samp, sigma = s))
  }

  return(sigma)
}

# ---------------------------------------------------------------------------
# Compute epsilon_j (cross-contamination rate) from Input/Nascent pairs.
# epsilon = (endo_reads_in_input) / (endo_reads_in_nascent + endo_reads_in_input)
# If no paired data is provided, epsilon defaults to 0.
# ---------------------------------------------------------------------------
compute_epsilon <- function(endo_counts, input_counts, nascent_counts) {
  endo_dt    <- endo_counts$data
  input_dt   <- input_counts$data
  nascent_dt <- nascent_counts$data
  gene_col   <- endo_counts$gene_col

  sample_cols <- setdiff(names(endo_dt), gene_col)

  epsilon <- data.table(sample = character(), epsilon = numeric())

  for (samp in sample_cols) {
    have_input   <- samp %in% names(input_dt)
    have_nascent <- samp %in% names(nascent_dt)

    if (!have_input || !have_nascent) {
      log_warn(sprintf("Sample '%s': missing input/nascent data; setting epsilon = 0", samp))
      epsilon <- rbind(epsilon, data.table(sample = samp, epsilon = 0))
      next
    }

    input_total   <- sum(input_dt[[samp]],   na.rm = TRUE)
    nascent_total <- sum(nascent_dt[[samp]], na.rm = TRUE)
    denom <- nascent_total + input_total

    if (denom == 0) {
      epsilon <- rbind(epsilon, data.table(sample = samp, epsilon = 0))
      next
    }

    e <- input_total / denom
    log_info(sprintf("  %s: epsilon = %.6f (input=%d, nascent=%d)", samp, e, input_total, nascent_total))
    epsilon <- rbind(epsilon, data.table(sample = samp, epsilon = e))
  }

  return(epsilon)
}

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
main <- function() {
  parser <- ArgumentParser(
    description = "Compute spike-in normalization factors for PRO-seq data"
  )
  parser$add_argument("--endogenous",
    required = TRUE,
    help = "Path to endogenous (gene body) count matrix (tab-separated)")
  parser$add_argument("--spike",
    required = TRUE,
    help = "Path to spike-in (dm6/external) count matrix (tab-separated)")
  parser$add_argument("--output",
    required = TRUE,
    help = "Path to output tab-separated factors file")
  parser$add_argument("--epsilon_input",
    default = NULL,
    help = "Path to input (No-4sU) count matrix for epsilon calculation (optional)")
  parser$add_argument("--epsilon_nascent",
    default = NULL,
    help = "Path to nascent (4sU-enriched) count matrix for epsilon calculation (optional)")

  args <- parser$parse_args()

  log_info("=== Spike Correction Factor Computation ===")

  # Read endogenous and spike matrices
  endo_counts  <- read_counts(args$endogenous, "endogenous")
  spike_counts <- read_counts(args$spike,      "spike-in")

  # Compute sigma
  sigma <- compute_sigma(endo_counts, spike_counts)

  # Compute epsilon if input/nascent data provided
  epsilon <- data.table(sample = sigma$sample, epsilon = 0)
  if (!is.null(args$epsilon_input) && !is.null(args$epsilon_nascent)) {
    log_info("Computing epsilon (cross-contamination rate) from Input/4sU paired data")
    input_counts   <- read_counts(args$epsilon_input,   "Input (No-4sU)")
    nascent_counts <- read_counts(args$epsilon_nascent, "Nascent (4sU-enriched)")
    epsilon <- compute_epsilon(endo_counts, input_counts, nascent_counts)
  } else {
    log_info("No epsilon data provided; all epsilon values set to 0")
  }

  # Merge and write output
  result <- merge(sigma, epsilon, by = "sample", all = TRUE)
  setnafill(result, fill = 0, cols = c("sigma", "epsilon"))

  log_info(sprintf("Writing output to: %s", args$output))
  fwrite(result, file = args$output, sep = "\t", quote = FALSE, row.names = FALSE)

  log_info("Done.")
}

# ---------------------------------------------------------------------------
if (sys.nframe() == 0) {
  main()
}
