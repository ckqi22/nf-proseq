#!/usr/bin/env Rscript
# =============================================================================
# pausing_index.R — Compute the PRO-seq pausing index for each gene
#
# PI_i = (promoter reads / promoter window length) / (genebody reads / genebody length)
#
# Input: two count matrices produced by the quantification subworkflow
#   (featurecounts_merge.R), each with columns: gene_id, length, <sample1>, ...
#   - promoter matrix:  'length' = promoter window size (constant across genes)
#   - genebody matrix:  'length' = genebody length (drives the min-length filter)
# Output: a single tab-separated PI table.
# =============================================================================
suppressWarnings(suppressMessages({
    library(argparser)
}))

# ── Parse args ──
argv <- arg_parser("Compute length-normalized PRO-seq pausing index")
argv <- add_argument(argv, "--promoter_count",  help = "Promoter count matrix (gene_id, length, samples)")
argv <- add_argument(argv, "--genebody_count",  help = "Genebody count matrix (gene_id, length, samples)")
argv <- add_argument(argv, "--min_gene_length", help = "Minimum gene body length (bp)", default = 800, type = "integer")
argv <- add_argument(argv, "--output",          help = "Output PI table (tsv)")
argv <- parse_args(argv)

# ── Compute pausing index ──
compute_pausing_index <- function(promoter, genebody, min_len) {
    gene_col <- colnames(promoter)[1]
    if (colnames(genebody)[1] != gene_col)
        stop("promoter and genebody count matrices must share the same first column (gene_id)")
    if (!"length" %in% colnames(promoter)) stop("promoter matrix must have a 'length' column")
    if (!"length" %in% colnames(genebody)) stop("genebody matrix must have a 'length' column")

    # Align both matrices on common genes, then filter by gene body length.
    common <- intersect(promoter[[gene_col]], genebody[[gene_col]])
    promoter <- promoter[match(common, promoter[[gene_col]]), , drop = FALSE]
    genebody <- genebody[match(common, genebody[[gene_col]]), , drop = FALSE]

    genebody_len <- as.numeric(genebody[["length"]])
    keep   <- genebody_len >= min_len
    promoter <- promoter[keep, , drop = FALSE]
    genebody <- genebody[keep, , drop = FALSE]

    genes   <- promoter[[gene_col]]
    promoter_len <- as.numeric(promoter[["length"]])   # promoter window length
    genebody_len <- as.numeric(genebody[["length"]])   # genebody length

    sample_cols <- setdiff(colnames(promoter), c(gene_col, "length"))

    out <- data.frame(gene_id = genes, promoter_len = promoter_len, genebody_len = genebody_len, stringsAsFactors = FALSE)
    for (s in sample_cols) {
        pause <- as.numeric(promoter[[s]])
        body  <- as.numeric(genebody[[s]])
        pi    <- (pause / promoter_len) / (body / genebody_len)
        pi[body <= 0] <- NA   # undefined when the gene body has no reads
        out[[paste0(s, "_promoter")]] <- pause
        out[[paste0(s, "_genebody")]]  <- body
        out[[paste0(s, "_PI")]]    <- pi
    }
    out
}

# ── Main ──
main <- function(argv) {
    if (!file.exists(argv$promoter_count)) stop("Promoter counts file not found: ", argv$promoter_count)
    if (!file.exists(argv$genebody_count)) stop("Genebody counts file not found: ", argv$genebody_count)

    promoter <- read.delim(argv$promoter_count, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
    genebody <- read.delim(argv$genebody_count, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)

    out <- compute_pausing_index(promoter, genebody, as.integer(argv$min_gene_length))

    write.table(out, file = argv$output, sep = "\t", quote = FALSE, row.names = FALSE)
    message("[pausing_index] ", nrow(out), " genes x ", ncol(out), " cols -> ", argv$output)
}

main(argv)
