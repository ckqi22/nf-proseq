#!/usr/bin/env Rscript
#
# metagene_matrix.R - Wrapper around deepTools computeMatrix, plotProfile, and plotHeatmap
#
# Runs strand-specific metagene analysis on PRO-seq plus/minus BAM files.
# Generates TSS and TES metagene PDFs with separate profiles for each strand.
#
# Requirements:
#   deepTools must be installed and available on PATH.
#
# Usage:
#   Rscript metagene_matrix.R --plus_bams bam1_plus.bam,bam2_plus.bam \
#       --minus_bams bam1_minus.bam,bam2_minus.bam \
#       --regions metagene_tss.bed --window 3000 --output_prefix my_meta --mode tss
#

suppressPackageStartupMessages({
  library(argparse)
  library(data.table)
  library(parallel)
})

# ---------------------------------------------------------------------------
# Logging
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
# Check deepTools availability
# ---------------------------------------------------------------------------
check_deeptools <- function() {
  ct_avail <- Sys.which("computeMatrix")
  pp_avail <- Sys.which("plotProfile")
  ph_avail <- Sys.which("plotHeatmap")

  if (ct_avail == "" || pp_avail == "") {
    log_error("deepTools (computeMatrix / plotProfile) not found on PATH.")
    log_error("Install with: conda install -c bioconda deeptools")
    quit(status = 1)
  }
  log_info(sprintf("deepTools found: computeMatrix=%s, plotProfile=%s, plotHeatmap=%s",
    ct_avail, pp_avail, ph_avail))
}

# ---------------------------------------------------------------------------
# Run deepTools computeMatrix for a set of BAM files
# ---------------------------------------------------------------------------
run_compute_matrix <- function(bam_files, regions, window_bp, output_prefix, mode, label) {
  bam_str <- paste(bam_files, collapse = " ")

  # Determine scale-regions vs reference-point based on mode
  if (mode == "tss") {
    matrix_cmd <- sprintf(
      'computeMatrix reference-point --referencePoint TSS -b %d -a %d -S %s -R %s -o %s_%s_matrix.gz --missingDataAsZero -p max 2>&1',
      window_bp, window_bp, bam_str, regions, output_prefix, label
    )
  } else if (mode == "tes") {
    matrix_cmd <- sprintf(
      'computeMatrix reference-point --referencePoint TES -b %d -a %d -S %s -R %s -o %s_%s_matrix.gz --missingDataAsZero -p max 2>&1',
      window_bp, window_bp, bam_str, regions, output_prefix, label
    )
  } else {
    log_error(sprintf("Unknown mode '%s'. Must be 'tss' or 'tes'.", mode))
    quit(status = 1)
  }

  log_info(sprintf("Running computeMatrix (%s):", label))
  log_info(sprintf("  CMD: %s", matrix_cmd))

  ret <- system(matrix_cmd, intern = TRUE)
  # Print output (useful for debugging)
  last_lines <- tail(ret, 5)
  for (l in last_lines) {
    cat(sprintf("  [computeMatrix] %s\n", l))
  }

  matrix_file <- sprintf("%s_%s_matrix.gz", output_prefix, label)
  if (!file.exists(matrix_file)) {
    log_error(sprintf("computeMatrix failed: %s not created.", matrix_file))
    quit(status = 1)
  }

  return(matrix_file)
}

# ---------------------------------------------------------------------------
# Run deepTools plotProfile
# ---------------------------------------------------------------------------
run_plot_profile <- function(matrix_file, output_prefix, label) {
  plot_file <- sprintf("%s_%s_profile.pdf", output_prefix, label)

  cmd <- sprintf(
    'plotProfile -m %s -o %s --perGroup --plotTitle "PRO-seq %s Metagene Profile - %s strand" --averageType mean --colors %s 2>&1',
    matrix_file, plot_file, toupper(substr(label, 1, 3)), label, label
  )

  # Use appropriate colors: blue for plus, red for minus
  if (label == "plus") {
    cmd <- sprintf(
      'plotProfile -m %s -o %s --perGroup --plotTitle "PRO-seq Metagene Profile - Plus Strand" --averageType mean 2>&1',
      matrix_file, plot_file
    )
  } else {
    cmd <- sprintf(
      'plotProfile -m %s -o %s --perGroup --plotTitle "PRO-seq Metagene Profile - Minus Strand" --averageType mean 2>&1',
      matrix_file, plot_file
    )
  }

  log_info(sprintf("Running plotProfile (%s): %s", label, cmd))

  ret <- system(cmd, intern = TRUE)
  for (l in tail(ret, 5)) {
    cat(sprintf("  [plotProfile] %s\n", l))
  }

  if (!file.exists(plot_file)) {
    log_error(sprintf("plotProfile failed: %s not created.", plot_file))
  } else {
    log_info(sprintf("Profile plot saved: %s", plot_file))
  }

  return(plot_file)
}

# ---------------------------------------------------------------------------
# Run deepTools plotHeatmap (optional)
# ---------------------------------------------------------------------------
run_plot_heatmap <- function(matrix_file, output_prefix, label, sort_regions = "descend") {
  heatmap_file <- sprintf("%s_%s_heatmap.pdf", output_prefix, label)

  cmd <- sprintf(
    'plotHeatmap -m %s -o %s --colorMap RdYlBu_r --sortRegions %s --whatToShow "plot, heatmap and colorbar" --missingDataColor white --plotTitle "PRO-seq Metagene Heatmap - %s Strand" 2>&1',
    matrix_file, heatmap_file, sort_regions, label
  )

  log_info(sprintf("Running plotHeatmap (%s): %s", label, cmd))

  ret <- system(cmd, intern = TRUE)
  for (l in tail(ret, 5)) {
    cat(sprintf("  [plotHeatmap] %s\n", l))
  }

  if (!file.exists(heatmap_file)) {
    log_error(sprintf("plotHeatmap failed: %s not created.", heatmap_file))
  } else {
    log_info(sprintf("Heatmap saved: %s", heatmap_file))
  }

  return(heatmap_file)
}

# ---------------------------------------------------------------------------
# Parse comma-separated list of BAM paths
# ---------------------------------------------------------------------------
parse_bam_list <- function(bam_string) {
  if (is.null(bam_string) || nchar(bam_string) == 0) {
    return(character(0))
  }
  bams <- unlist(strsplit(bam_string, ","))
  bams <- trimws(bams)

  # Verify existence
  missing <- bams[!file.exists(bams)]
  if (length(missing) > 0) {
    log_error(sprintf("BAM files not found: %s", paste(missing, collapse = ", ")))
    quit(status = 1)
  }
  log_info(sprintf("Parsed %d BAM files from: %s", length(bams), bam_string))
  return(bams)
}

# ---------------------------------------------------------------------------
# Validate BED file and check strand separation
# ---------------------------------------------------------------------------
validate_bed <- function(bed_path) {
  if (!file.exists(bed_path)) {
    log_error(sprintf("BED file not found: %s", bed_path))
    quit(status = 1)
  }
  n_lines <- length(readLines(bed_path))
  log_info(sprintf("BED file: %s (%d regions)", bed_path, n_lines))
  return(bed_path)
}

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
main <- function() {
  parser <- ArgumentParser(
    description = "deepTools wrapper for strand-specific PRO-seq metagene analysis"
  )
  parser$add_argument("--plus_bams",
    required = TRUE,
    help = "Comma-separated list of plus-strand BAM files")
  parser$add_argument("--minus_bams",
    required = TRUE,
    help = "Comma-separated list of minus-strand BAM files")
  parser$add_argument("--regions",
    required = TRUE,
    help = "BED file of regions (TSS or TES) for metagene analysis")
  parser$add_argument("--window",
    type    = "integer",
    default = 3000L,
    help    = "Window size in bp around reference point (default: 3000)")
  parser$add_argument("--output_prefix",
    required = TRUE,
    help = "Prefix for output files")
  parser$add_argument("--mode",
    required = TRUE,
    help = "Analysis mode: 'tss' or 'tes'")
  parser$add_argument("--heatmap",
    action = "store_true",
    default = FALSE,
    help = "Also generate heatmaps (default: profile only)")

  args <- parser$parse_args()

  args$mode <- tolower(args$mode)
  if (!args$mode %in% c("tss", "tes")) {
    log_error("--mode must be 'tss' or 'tes'")
    quit(status = 1)
  }

  log_info("=== PRO-seq Metagene Analysis ===")
  log_info(sprintf("Mode: %s, Window: +/-%d bp", args$mode, args$window))

  # Check dependencies
  check_deeptools()

  # Parse BAM lists
  plus_bams  <- parse_bam_list(args$plus_bams)
  minus_bams <- parse_bam_list(args$minus_bams)

  # Validate BED
  validate_bed(args$regions)

  # Run analysis for each strand
  outputs <- list()

  for (strand_label in c("plus", "minus")) {
    bam_files <- if (strand_label == "plus") plus_bams else minus_bams

    if (length(bam_files) == 0) {
      log_warn(sprintf("No %s-strand BAM files provided; skipping", strand_label))
      next
    }

    log_info(sprintf("--- %s strand analysis (%d BAMs) ---", strand_label, length(bam_files)))

    # computeMatrix
    matrix_file <- run_compute_matrix(
      bam_files      = bam_files,
      regions        = args$regions,
      window_bp      = args$window,
      output_prefix  = args$output_prefix,
      mode           = args$mode,
      label          = strand_label
    )

    # plotProfile
    profile_file <- run_plot_profile(matrix_file, args$output_prefix, strand_label)
    outputs[[paste0(strand_label, "_profile")]] <- profile_file

    # plotHeatmap (optional)
    if (args$heatmap) {
      heatmap_file <- run_plot_heatmap(matrix_file, args$output_prefix, strand_label)
      outputs[[paste0(strand_label, "_heatmap")]] <- heatmap_file
    }
  }

  log_info("--- Summary ---")
  for (name in names(outputs)) {
    log_info(sprintf("  %s: %s", name, outputs[[name]]))
  }
  log_info("Done.")
}

# ---------------------------------------------------------------------------
if (sys.nframe() == 0) {
  main()
}
