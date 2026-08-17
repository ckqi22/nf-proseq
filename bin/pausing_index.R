#!/usr/bin/env Rscript
# =============================================================================
# pausing_index.R — Compute the PRO-seq pausing index for each gene
#
# PI_i = (TSS reads / TSS window length) / (gene body reads / gene body length)
#
# Input: two count matrices produced by the quantification subworkflow
#   (featurecounts_merge.R), each with columns: gene_id, length, <sample1>, ...
#   - TSS matrix:       'length' = TSS window size (constant across genes)
#   - gene body matrix: 'length' = gene body length (drives the min-length filter)
# Output: a single tab-separated PI table.
# =============================================================================
suppressWarnings(suppressMessages({
    library(argparser)
}))

# ── Parse args ──
argv <- arg_parser("Compute length-normalized PRO-seq pausing index")
argv <- add_argument(argv, "--tss_counts",      help = "TSS count matrix (gene_id, length, samples)")
argv <- add_argument(argv, "--gb_counts",       help = "Gene body count matrix (gene_id, length, samples)")
argv <- add_argument(argv, "--min_gene_length", help = "Minimum gene body length (bp)", default = 800, type = "integer")
argv <- add_argument(argv, "--output",          help = "Output PI table (tsv)")
argv <- parse_args(argv)

# ── Compute pausing index ──
compute_pausing_index <- function(tss, gb, min_len) {
    gene_col <- colnames(tss)[1]
    if (colnames(gb)[1] != gene_col)
        stop("TSS and gene body matrices must share the same first column (gene_id)")
    if (!"length" %in% colnames(tss)) stop("TSS matrix must have a 'length' column")
    if (!"length" %in% colnames(gb))  stop("gene body matrix must have a 'length' column")

    # Align both matrices on common genes, then filter by gene body length.
    common <- intersect(tss[[gene_col]], gb[[gene_col]])
    tss <- tss[match(common, tss[[gene_col]]), , drop = FALSE]
    gb  <- gb[match(common, gb[[gene_col]]), , drop = FALSE]

    gb_len <- as.numeric(gb[["length"]])
    keep   <- gb_len >= min_len
    tss <- tss[keep, , drop = FALSE]
    gb  <- gb[keep, , drop = FALSE]

    genes   <- tss[[gene_col]]
    tss_len <- as.numeric(tss[["length"]])   # TSS window length
    gb_len  <- as.numeric(gb[["length"]])    # gene body length

    sample_cols <- setdiff(colnames(tss), c(gene_col, "length"))

    out <- data.frame(gene_id = genes, tss_len = tss_len, genebody_len = gb_len, stringsAsFactors = FALSE)
    for (s in sample_cols) {
        pause <- as.numeric(tss[[s]])
        body  <- as.numeric(gb[[s]])
        pi    <- (pause / tss_len) / (body / gb_len)
        pi[body <= 0] <- NA   # undefined when the gene body has no reads
        out[[paste0(s, "_tss")]] <- pause
        out[[paste0(s, "_genebody")]]  <- body
        out[[paste0(s, "_PI")]]    <- pi
    }
    out
}

# ── Main ──
main <- function(argv) {
    if (!file.exists(argv$tss_counts)) stop("TSS counts file not found: ", argv$tss_counts)
    if (!file.exists(argv$gb_counts))  stop("gene body counts file not found: ", argv$gb_counts)

    tss <- read.delim(argv$tss_counts, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
    gb  <- read.delim(argv$gb_counts,  header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)

    out <- compute_pausing_index(tss, gb, as.integer(argv$min_gene_length))

    write.table(out, file = argv$output, sep = "\t", quote = FALSE, row.names = FALSE)
    message("[pausing_index] ", nrow(out), " genes x ", ncol(out), " cols -> ", argv$output)
}

main(argv)
