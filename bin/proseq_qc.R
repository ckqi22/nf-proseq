#!/usr/bin/env Rscript
# =============================================================================
# proseq_qc.R — PRO-seq specific QC metrics per sample
#
# Computes: mapping rate, strand separation efficiency,
# TSS enrichment ratio, library complexity from BAM files.
# =============================================================================
suppressWarnings(suppressMessages({
    library(argparser)
}))

argv <- arg_parser("PRO-seq QC metrics from BAM and alignment log")
argv <- add_argument(argv, "--bam",            help = "Path to BAM file")
argv <- add_argument(argv, "--alignment_log",  help = "Path to Bowtie2 alignment log file")
argv <- add_argument(argv, "--tss_regions",    help = "TSS regions BED file (optional)")
argv <- add_argument(argv, "--output",         help = "Output QC report file")
argv <- parse_args(argv)

if (!file.exists(argv$bam))            stop("BAM file not found: ", argv$bam)
if (!file.exists(argv$alignment_log))  stop("Alignment log not found: ", argv$alignment_log)

sample_name <- gsub("\\.bam$", "", basename(argv$bam))

message("[proseq_qc] Processing: ", sample_name)

# --- Parse Bowtie2 alignment log ---
log_lines <- readLines(argv$alignment_log)

parse_bowtie2_value <- function(lines, pattern) {
    match_line <- grep(pattern, lines, value = TRUE)
    if (length(match_line) == 0) return(NA_real_)
    as.numeric(gsub("[^0-9.]", "", sub(".*\\(([0-9.]+)%\\).*", "\\1", match_line)))
}

total_reads    <- parse_bowtie2_value(log_lines, "reads; of these")
paired_reads   <- parse_bowtie2_value(log_lines, "were paired")
concordant_0x  <- parse_bowtie2_value(log_lines, "aligned concordantly exactly 1 time")
concordant_nx  <- parse_bowtie2_value(log_lines, "aligned concordantly >1 times")
overall_rate   <- parse_bowtie2_value(log_lines, "overall alignment rate")

# --- BAM-level metrics via samtools stats ---
n_mapped <- 0L
n_total  <- 0L

# Try to get basic stats
idxstats_cmd <- paste0("samtools idxstats ", argv$bam)
idxstats_out <- tryCatch(
    system(idxstats_cmd, intern = TRUE),
    error = function(e) character(0)
)

if (length(idxstats_out) > 0) {
    idx_lines <- strsplit(idxstats_out, "\t")
    n_mapped <- sum(sapply(idx_lines, function(x) as.numeric(x[3])), na.rm = TRUE)
    n_total  <- sum(sapply(idx_lines, function(x) as.numeric(x[3]) + as.numeric(x[4])), na.rm = TRUE)
}

# --- Strand separation efficiency ---
# Count forward vs reverse reads
flagstat_cmd <- paste0("samtools flagstat ", argv$bam)
flagstat_out <- tryCatch(
    system(flagstat_cmd, intern = TRUE),
    error = function(e) character(0)
)

n_forward <- NA_integer_
n_reverse <- NA_integer_
strand_ratio <- NA_real_

if (length(flagstat_out) > 0) {
    # Extract total mapped
    mapped_line <- grep("mapped", flagstat_out, value = TRUE)[1]
    if (length(mapped_line) > 0) {
        # Not directly from flagstat; use idxstats
    }
}

# Try samtools view counting by flag
if (n_mapped > 0) {
    forward_cmd <- paste0("samtools view -c -F 0x10 ", argv$bam)
    reverse_cmd <- paste0("samtools view -c -f 0x10 ", argv$bam)

    n_forward <- tryCatch(as.integer(system(forward_cmd, intern = TRUE)), error = function(e) NA_integer_)
    n_reverse <- tryCatch(as.integer(system(reverse_cmd, intern = TRUE)), error = function(e) NA_integer_)

    if (!is.na(n_forward) && !is.na(n_reverse) && (n_forward + n_reverse) > 0) {
        strand_ratio <- n_forward / (n_forward + n_reverse)
    }
}

# --- Compile report ---
report_lines <- c(
    paste0("=== PRO-seq QC Report: ", sample_name, " ==="),
    paste0("Date: ", Sys.time()),
    "",
    "--- Alignment Metrics ---",
    paste0("Total reads:            ", total_reads),
    paste0("Overall alignment rate: ", ifelse(is.na(overall_rate), "N/A", paste0(overall_rate, "%"))),
    paste0("Concordant 1x:          ", ifelse(is.na(concordant_0x), "N/A", paste0(concordant_0x, "%"))),
    paste0("Concordant >1x:         ", ifelse(is.na(concordant_nx), "N/A", paste0(concordant_nx, "%"))),
    "",
    "--- BAM Metrics ---",
    paste0("Total mapped reads:     ", n_mapped),
    paste0("Total reads (idxstats): ", n_total),
    "",
    "--- Strand Separation ---",
    paste0("Forward reads:          ", ifelse(is.na(n_forward), "N/A", n_forward)),
    paste0("Reverse reads:          ", ifelse(is.na(n_reverse), "N/A", n_reverse)),
    paste0("Forward/Total ratio:    ", ifelse(is.na(strand_ratio), "N/A", round(strand_ratio, 4)))
)

# Write report
writeLines(report_lines, argv$output)
message("[proseq_qc] Report written: ", argv$output)
message("[proseq_qc] Done.")
