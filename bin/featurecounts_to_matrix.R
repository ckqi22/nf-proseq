#!/usr/bin/env Rscript
# =============================================================================
# featurecounts_to_matrix.R — Convert a featureCounts count file into a clean
# count matrix for pausing-index computation.
#
# featureCounts output:  Geneid  Chr  Start  End  Strand  Length  <sample>.bam ...
# Clean matrix output:   gene_id  length  <sample1>  <sample2> ...
#
#   - skip the leading '# Program:featureCounts ...' comment line(s)
#   - drop non-count columns (Chr / Start / End / Strand)
#   - rename Geneid -> gene_id, Length -> length (gene body length)
#   - strip trailing '.bam' from sample columns so they match group/sample names
# =============================================================================
suppressWarnings(suppressMessages({
    library(argparser)
}))

# ── Parse args ──
argv <- arg_parser("Convert a featureCounts count file into a clean count matrix")
argv <- add_argument(argv, "--input",  help = "featureCounts count file")
argv <- add_argument(argv, "--output", help = "Output clean count matrix (tsv)")
argv <- parse_args(argv)

# ── Main ──
main <- function(argv) {
    # comment.char="#" skips the featureCounts '# Program:...' header line(s)
    dt <- read.delim(argv$input, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
                     comment.char = "#", check.names = FALSE)

    names(dt)[names(dt) == "Geneid"] <- "gene_id"
    names(dt)[names(dt) == "Length"] <- "length"

    drop <- intersect(c("Chr", "Start", "End", "Strand"), names(dt))
    if (length(drop) > 0) {
        dt <- dt[, setdiff(names(dt), drop), drop = FALSE]
    }

    sample_cols <- setdiff(names(dt), c("gene_id", "length"))
    if (length(sample_cols) > 0) {
        names(dt)[match(sample_cols, names(dt))] <- sub("\\.bam$", "", sample_cols, ignore.case = TRUE)
    }

    write.table(dt, file = argv$output, sep = "\t", quote = FALSE, row.names = FALSE)
    message("[featurecounts_to_matrix] ", nrow(dt), " genes x ", ncol(dt), " cols -> ", argv$output)
}

main(argv)
