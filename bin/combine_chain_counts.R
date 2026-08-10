#!/usr/bin/env Rscript
# =============================================================================
# combine_chain_counts.R — Merge plus/minus strand featureCounts outputs
#
# For PRO-seq strand-specific analysis:
#   - Plus-strand BAM reads → minus-strand gene counts
#   - Minus-strand BAM reads → plus-strand gene counts
# This script combines the two count files into a single gene-level matrix.
# =============================================================================
suppressWarnings(suppressMessages({
    library(argparser)
}))

argv <- arg_parser("Merge plus-strand and minus-strand featureCounts outputs")
argv <- add_argument(argv, "--plus",   help = "featureCounts output for plus-strand BAM")
argv <- add_argument(argv, "--minus",  help = "featureCounts output for minus-strand BAM")
argv <- add_argument(argv, "--output", help = "Output combined counts file (tab-separated)")
argv <- parse_args(argv)

# Validate inputs
if (!file.exists(argv$plus))  stop("Plus-strand counts file not found: ", argv$plus)
if (!file.exists(argv$minus)) stop("Minus-strand counts file not found: ", argv$minus)

message("[combine_chain_counts] Reading plus-strand counts: ", argv$plus)
plus <- read.delim(argv$plus, header = TRUE, stringsAsFactors = FALSE, comment.char = "#")
message("  ", nrow(plus), " genes x ", ncol(plus), " columns")

message("[combine_chain_counts] Reading minus-strand counts: ", argv$minus)
minus <- read.delim(argv$minus, header = TRUE, stringsAsFactors = FALSE, comment.char = "#")
message("  ", nrow(minus), " genes x ", ncol(minus), " columns")

# featureCounts output columns: Geneid, Chr, Start, End, Strand, Length, <sample>.bam
# The count column is typically column 7
count_col_plus  <- ncol(plus)
count_col_minus <- ncol(minus)

# Build combined table: take annotation from plus file, add both counts
combined <- data.frame(
    Geneid        = plus$Geneid,
    Chr           = plus$Chr,
    Start         = plus$Start,
    End           = plus$End,
    Strand        = plus$Strand,
    Length        = plus$Length,
    Plus_Counts   = plus[, count_col_plus],
    Minus_Counts  = 0,
    Total_Counts  = 0,
    stringsAsFactors = FALSE
)

# Match minus counts by Geneid
minus_map <- setNames(minus[, count_col_minus], minus$Geneid)
matched_idx <- match(combined$Geneid, names(minus_map))
combined$Minus_Counts <- ifelse(is.na(matched_idx), 0, minus_map[matched_idx])
combined$Total_Counts  <- combined$Plus_Counts + combined$Minus_Counts

# Rename count columns to use sample name from the input filename
plus_sample  <- gsub("\\.gene_body_plus\\.counts\\.txt$", "", basename(argv$plus))
plus_sample  <- gsub("\\.tss_plus\\.counts\\.txt$", "", plus_sample)
minus_sample <- gsub("\\.gene_body_minus\\.counts\\.txt$", "", basename(argv$minus))
minus_sample <- gsub("\\.tss_minus\\.counts\\.txt$", "", minus_sample)

names(combined)[names(combined) == "Plus_Counts"]  <- paste0(plus_sample, "_Plus")
names(combined)[names(combined) == "Minus_Counts"] <- paste0(minus_sample, "_Minus")
names(combined)[names(combined) == "Total_Counts"] <- paste0(plus_sample, "_Total")

message("[combine_chain_counts] Writing output: ", argv$output)
write.table(combined, file = argv$output, sep = "\t", quote = FALSE, row.names = FALSE)
message("[combine_chain_counts] Done. ", nrow(combined), " genes written.")
