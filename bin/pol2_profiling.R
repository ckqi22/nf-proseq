#!/usr/bin/env Rscript
#
# pol2_profiling.R - Single-nucleotide resolution Pol II occupancy profiling from PRO-seq BAM files
#
# Extracts 5' end positions of reads from BAM files, counts them at single-nucleotide
# resolution genome-wide, normalizes by spike-in factors or total library size,
# and produces per-sample signal tracks.
#
# Usage:
#   Rscript pol2_profiling.R --bam_dir <dir> --spike_factors <factors.txt> --output_dir <dir>
#

suppressPackageStartupMessages({
  library(argparse)
  library(data.table)
  library(GenomicRanges)
  library(Rsamtools)
  library(GenomicAlignments)
  library(rtracklayer)
  library(openxlsx)
  library(ComplexHeatmap)
  library(circlize)
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
# Discover BAM files in a directory
# ---------------------------------------------------------------------------
discover_bam_files <- function(bam_dir) {
  if (!dir.exists(bam_dir)) {
    log_error(sprintf("BAM directory not found: %s", bam_dir))
    quit(status = 1)
  }
  bam_files <- list.files(bam_dir, pattern = "\\.bam$", full.names = TRUE, ignore.case = TRUE)
  if (length(bam_files) == 0) {
    log_error(sprintf("No BAM files found in: %s", bam_dir))
    quit(status = 1)
  }
  log_info(sprintf("Found %d BAM files:", length(bam_files)))
  for (bf in bam_files) {
    log_info(sprintf("  %s", basename(bf)))
  }
  return(bam_files)
}

# ---------------------------------------------------------------------------
# Read spike-in factors
# ---------------------------------------------------------------------------
read_spike_factors <- function(path) {
  if (is.null(path) || !file.exists(path)) {
    log_warn("No spike factors file; using library-size normalization")
    return(NULL)
  }
  sf <- fread(path, header = TRUE, sep = "\t")
  log_info(sprintf("Read spike factors for %d samples", nrow(sf)))
  return(sf)
}

# ---------------------------------------------------------------------------
# Extract 5' ends from a BAM file and count per position
# For single-end PRO-seq: the 5' end of the read = the 5' end of the nascent RNA.
# For paired-end: use the 5' end of read1 (R1) for plus strand, read2 (R2) for minus strand.
# This implementation handles single-end reads.
# ---------------------------------------------------------------------------
extract_five_prime_counts <- function(bam_path, sample_name) {
  log_info(sprintf("  Processing: %s", basename(bam_path)))

  # Read BAM as GAlignments
  param <- ScanBamParam(
    what = c("rname", "pos", "strand", "qwidth"),
    flag  = scanBamFlag(isUnmappedQuery = FALSE)
  )

  ga <- readGAlignments(bam_path, param = param)

  if (length(ga) == 0) {
    log_warn(sprintf("  No reads in %s", basename(bam_path)))
    return(data.table(chr = character(), pos = integer(), strand = character(), count = integer()))
  }

  # Determine 5' end position based on strand
  # For plus strand reads: 5' end = start(ga)
  # For minus strand reads: 5' end = end(ga)
  strand_vec <- as.character(strand(ga))

  five_prime_pos <- ifelse(strand_vec == "+",
    start(ga),   # plus strand: 5' is at start
    end(ga)      # minus strand: 5' is at end
  )

  dt <- data.table(
    chr     = as.character(seqnames(ga)),
    pos     = five_prime_pos,
    strand  = strand_vec
  )

  # Count occurrences of each position
  counts <- dt[, .(count = .N), by = .(chr, pos, strand)]
  setnames(counts, c("chr", "pos", "strand", sample_name))

  log_info(sprintf("    %d unique 5' positions from %d reads", nrow(counts), length(ga)))

  return(counts)
}

# ---------------------------------------------------------------------------
# Normalize counts by spike-in factor or library size
# Normalized signal = count / sigma_j  (if spike factors available)
#                   = count / (total_counts / 1e6)  (CPM otherwise)
# ---------------------------------------------------------------------------
normalize_counts <- function(counts_list, spike_factors) {
  all_samples <- names(counts_list)
  norm_list <- list()

  for (samp in all_samples) {
    dt <- counts_list[[samp]]
    raw_total <- sum(dt[[samp]], na.rm = TRUE)

    if (!is.null(spike_factors) && samp %in% spike_factors$sample) {
      sigma <- spike_factors[sample == samp, ]$sigma
      if (!is.na(sigma) && sigma > 0) {
        norm_factor <- sigma
        label <- "spike"
      } else {
        norm_factor <- raw_total / 1e6
        label <- "CPM"
      }
    } else {
      norm_factor <- raw_total / 1e6
      label <- "CPM"
    }

    norm_col <- paste0(samp, "_norm")
    dt[[norm_col]] <- dt[[samp]] / norm_factor

    norm_list[[samp]] <- dt
    log_info(sprintf("  %s: total=%d, norm=%s, factor=%.6f", samp, raw_total, label, norm_factor))
  }

  return(norm_list)
}

# ---------------------------------------------------------------------------
# Merge per-sample counts into a single genome-wide table
# ---------------------------------------------------------------------------
merge_genome_wide <- function(counts_list, spike_factors) {
  norm_list <- normalize_counts(counts_list, spike_factors)

  # Full outer join all samples by chr, pos, strand
  result <- NULL
  for (samp in names(norm_list)) {
    dt <- norm_list[[samp]]
    if (is.null(result)) {
      result <- dt
    } else {
      result <- merge(result, dt, by = c("chr", "pos", "strand"), all = TRUE)
    }
  }

  # Replace NAs with 0
  for (col in names(result)) {
    if (col %in% c("chr", "pos", "strand")) next
    setnafill(result, fill = 0, cols = col)
  }

  # Add start/end columns (1-based single nucleotide)
  result[, start := pos]
  result[, end   := pos]

  # Reorder columns
  col_order <- c("chr", "start", "end", "strand",
    intersect(names(result), c("chr", "start", "end", "strand", "pos")),
    setdiff(names(result), c("chr", "start", "end", "strand", "pos")))

  # Sort by chr then pos
  setorder(result, chr, pos)

  log_info(sprintf("Merged genome-wide table: %d positions", nrow(result)))
  return(result)
}

# ---------------------------------------------------------------------------
# Generate a heatmap of TSS-region signal (requires a TSS BED file)
# ---------------------------------------------------------------------------
generate_tss_heatmap <- function(counts_list, tss_bed_path, output_dir, window = 3000) {
  if (!file.exists(tss_bed_path)) {
    log_warn(sprintf("TSS BED file not found: %s; skipping heatmap", tss_bed_path))
    return(NULL)
  }

  tss_regions <- fread(tss_bed_path, header = FALSE, sep = "\t")
  if (ncol(tss_regions) < 6) {
    log_warn("TSS BED must have at least 6 columns (chrom, start, end, name, score, strand); skipping")
    return(NULL)
  }
  setnames(tss_regions, c("chr", "start", "end", "name", "score", "strand"))

  log_info(sprintf("Generating TSS heatmap for %d genes...", nrow(tss_regions)))

  # For each TSS, extract signal in a surrounding window
  n_genes <- min(nrow(tss_regions), 500)  # cap for visualization
  half_win <- floor(window / 2)
  mat_list <- list()

  for (samp in names(counts_list)) {
    dt <- counts_list[[samp]]
    mat <- matrix(NA, nrow = n_genes, ncol = window)

    for (i in seq_len(n_genes)) {
      region <- tss_regions[i, ]
      tss_pos <- ifelse(region$strand == "+", region$start, region$end)
      region_chr <- region$chr
      win_start <- tss_pos - half_win
      win_end   <- tss_pos + half_win - 1

      # Extract counts in this window
      sub <- dt[chr == region_chr & pos >= win_start & pos <= win_end, ]
      if (nrow(sub) > 0) {
        idx <- sub$pos - win_start + 1
        # Ensure index is within bounds
        valid_idx <- idx >= 1 & idx <= window
        mat[i, idx[valid_idx]] <- sub[[samp]][valid_idx]
      }
      mat[i, is.na(mat[i, ])] <- 0
    }
    mat_list[[samp]] <- mat
  }

  # Average across samples or create multi-panel heatmap
  # For simplicity, plot the first sample
  samp_name <- names(counts_list)[1]
  mat <- mat_list[[samp_name]]

  # Log-transform for visualization
  mat_log <- log2(mat + 1)

  col_fun <- colorRamp2(
    breaks = c(0, quantile(mat_log[mat_log > 0], 0.5, na.rm = TRUE),
               quantile(mat_log, 0.95, na.rm = TRUE)),
    colors = c("white", "yellow", "red3")
  )

  pdf_path <- file.path(output_dir, "PROSeq_TSS_Heatmap.pdf")
  pdf(pdf_path, width = 10, height = 8)

  ht <- Heatmap(mat_log,
    name                = sprintf("log2(counts+1)\n%s", samp_name),
    col                 = col_fun,
    cluster_rows        = TRUE,
    cluster_columns     = FALSE,
    show_row_names      = FALSE,
    show_column_names   = FALSE,
    column_title        = sprintf("Distance from TSS (bp) - %s", samp_name),
    row_title           = "Genes",
    use_raster           = TRUE,
    heatmap_legend_param = list(direction = "horizontal")
  )
  draw(ht, heatmap_legend_side = "bottom")

  dev.off()
  log_info(sprintf("TSS heatmap saved: %s", pdf_path))
  return(pdf_path)
}

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
main <- function() {
  parser <- ArgumentParser(
    description = "Single-nucleotide Pol II occupancy profiling from PRO-seq BAM files"
  )
  parser$add_argument("--bam_dir",
    required = TRUE,
    help = "Directory containing BAM files")
  parser$add_argument("--spike_factors",
    default = NULL,
    help = "Spike correction factors file (tab-separated)")
  parser$add_argument("--output_dir",
    required = TRUE,
    help = "Output directory")
  parser$add_argument("--tss_bed",
    default = NULL,
    help = "TSS BED file for heatmap generation")

  args <- parser$parse_args()

  dir.create(args$output_dir, showWarnings = FALSE, recursive = TRUE)

  log_info("=== Pol II Single-Nucleotide Profiling ===")

  # Discover BAM files
  bam_files <- discover_bam_files(args$bam_dir)
  spike     <- read_spike_factors(args$spike_factors)

  # Extract 5' end counts per sample
  counts_list <- list()
  for (bam_path in bam_files) {
    sample_name <- sub("\\.bam$", "", basename(bam_path), ignore.case = TRUE)
    log_info(sprintf("Extracting 5' ends: %s", sample_name))
    counts_list[[sample_name]] <- extract_five_prime_counts(bam_path, sample_name)
  }

  # Merge and normalize
  genome_wide <- merge_genome_wide(counts_list, spike)

  # Write output
  xlsx_path <- file.path(args$output_dir, "PROSeq_profiling.xlsx")
  log_info(sprintf("Writting profiling data to: %s", xlsx_path))
  write.xlsx(genome_wide, file = xlsx_path, sheetName = "Pol2_Profiling", rowNames = FALSE)

  # Optional: write per-chromosome CSV files for large genomes
  bedgraph_dir <- file.path(args$output_dir, "bedgraph")
  dir.create(bedgraph_dir, showWarnings = FALSE)
  for (samp in names(counts_list)) {
    # Write raw counts as bedGraph
    norm_col <- paste0(samp, "_norm")
    cols_present <- intersect(c(samp, norm_col), names(genome_wide))
    for (col in cols_present) {
      bg_file <- file.path(bedgraph_dir, paste0(samp, "_", col, ".bedgraph"))
      bg_dt <- genome_wide[, .(chr, start, end, get(col))]
      fwrite(bg_dt, file = bg_file, sep = "\t", col.names = FALSE, quote = FALSE)
    }
  }
  log_info(sprintf("bedGraph files written to: %s", bedgraph_dir))

  # TSS heatmap if requested
  if (!is.null(args$tss_bed)) {
    generate_tss_heatmap(counts_list, args$tss_bed, args$output_dir)
  }

  log_info("Done.")
}

# ---------------------------------------------------------------------------
if (sys.nframe() == 0) {
  main()
}
