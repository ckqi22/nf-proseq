#!/usr/bin/env Rscript
# =============================================================================
# igv_tracks.R — Organize BAM and BigWig files for IGV visualization
#
# Creates organized directory with BAM/bigWig files and IGV batch script.
# =============================================================================
suppressWarnings(suppressMessages({
    library(argparser)
}))

argv <- arg_parser("Organize PRO-seq tracks for IGV visualization")
argv <- add_argument(argv, "--bam_dir",   help = "Directory containing BAM files")
argv <- add_argument(argv, "--bw_dir",    help = "Directory containing BigWig files")
argv <- add_argument(argv, "--output_dir", help = "Output directory for IGV tracks")
argv <- parse_args(argv)

dir.create(argv$output_dir, recursive = TRUE, showWarnings = FALSE)

message("[igv_tracks] Organizing IGV tracks...")

# --- Copy BAM files ---
bam_dest <- file.path(argv$output_dir, "bam")
dir.create(bam_dest, showWarnings = FALSE)

bam_files <- character(0)
if (dir.exists(argv$bam_dir)) {
    bam_files <- list.files(argv$bam_dir, pattern = "\\.bam$", full.names = TRUE, ignore.case = TRUE)
    for (bf in bam_files) {
        message("  Copying BAM: ", basename(bf))
        file.copy(bf, file.path(bam_dest, basename(bf)), overwrite = TRUE)
        # Also copy index
        bai_file <- paste0(bf, ".bai")
        if (file.exists(bai_file)) {
            file.copy(bai_file, file.path(bam_dest, basename(bai_file)), overwrite = TRUE)
        }
    }
}
message("[igv_tracks] ", length(bam_files), " BAM files")

# --- Copy BigWig files ---
bw_dest <- file.path(argv$output_dir, "bigwig")
dir.create(bw_dest, showWarnings = FALSE)

bw_files <- character(0)
if (dir.exists(argv$bw_dir)) {
    bw_files <- list.files(argv$bw_dir, pattern = "\\.(bw|bigwig|bigWig)$", full.names = TRUE)
    for (bf in bw_files) {
        message("  Copying BigWig: ", basename(bf))
        file.copy(bf, file.path(bw_dest, basename(bf)), overwrite = TRUE)
    }
}
message("[igv_tracks] ", length(bw_files), " BigWig files")

# --- Generate README ---
readme_lines <- c(
    "IGV Track Directory — PRO-seq",
    "==============================",
    "",
    "Directory structure:",
    "  bam/     — Strand-separated BAM files (+.bam = plus strand, -.bam = minus strand)",
    "  bigwig/  — Strand-separated BigWig coverage tracks (+.bw, -.bw)",
    "",
    "Loading in IGV:",
    "  1. Open IGV",
    "  2. Select the appropriate genome (e.g., hg38)",
    "  3. File > Load from File... > select .bw files from bigwig/",
    "  4. Color plus strand tracks blue and minus strand tracks red",
    "",
    paste0("Generated: ", Sys.time()),
    "",
    paste0("BAM files:  ", length(bam_files)),
    paste0("BigWig files: ", length(bw_files))
)
writeLines(readme_lines, file.path(argv$output_dir, "README.txt"))

message("[igv_tracks] Done: ", argv$output_dir)
