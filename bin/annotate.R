#!/usr/bin/env Rscript
# =============================================================================
# annotate.R — Merge raw counts + normalized values + gene annotation into one
# annotated profile table. Inner join on gene_id: only genes present in the
# annotation table are kept.
#
# Inputs:
#   --raw         raw count matrix       : gene_id  length  <sample>...
#   --normalized  normalize.R output     : gene_id  <s>.Fpkm ...  <s>.Cpm ...
#   --annotation  local annotation table : gene_id  gene_name ... GO_BP
#
# Output columns:
#   Gene_id  <s>.Count ...  <s>.Fpkm ...  <s>.Cpm ...  gene_name ... GO_BP
# =============================================================================
suppressWarnings(suppressMessages({
    library(argparser)
}))

# ── Parse args ──
argv <- arg_parser("Merge raw + normalized counts with gene annotation")
argv <- add_argument(argv, "--raw",        help = "Raw count matrix (gene_id, length, samples)")
argv <- add_argument(argv, "--normalized", help = "Normalized matrix (gene_id, <s>.<Suffix>)")
argv <- add_argument(argv, "--annotation", help = "Gene annotation table (first column = gene_id)")
argv <- add_argument(argv, "--output",     help = "Output annotated table (tsv)")
argv <- parse_args(argv)

main <- function(argv) {
    raw  <- read.delim(argv$raw,        header = TRUE, sep = "\t",
                       stringsAsFactors = FALSE, check.names = FALSE)
    norm <- read.delim(argv$normalized, header = TRUE, sep = "\t",
                       stringsAsFactors = FALSE, check.names = FALSE)
    ann  <- read.delim(argv$annotation, header = TRUE, sep = "\t",
                       stringsAsFactors = FALSE, check.names = FALSE,
                       quote = "", comment.char = "")
    ann[is.na(ann)] <- ""

    # raw -> Gene_id + <s>.Count (length dropped)
    sample_cols <- setdiff(names(raw), c("gene_id", "length"))
    if (length(sample_cols) == 0)
        stop("[annotate] no sample columns found in ", argv$raw)
    raw_out <- raw[, c("gene_id", sample_cols), drop = FALSE]
    names(raw_out) <- c("Gene_id", paste0(sample_cols, ".Count"))

    # normalized -> Gene_id + <s>.<Suffix>
    names(norm)[names(norm) == "gene_id"] <- "Gene_id"
    norm_cols <- setdiff(names(norm), "Gene_id")

    # annotation -> Gene_id + annotation columns (original order preserved)
    if (ncol(ann) < 2)
        stop("[annotate] annotation table has no annotation columns: ", argv$annotation)
    ann_cols <- names(ann)[-1]
    names(ann)[1] <- "Gene_id"

    # inner join on Gene_id: only genes present in the annotation table
    res <- merge(raw_out, norm, by = "Gene_id", all = FALSE, sort = FALSE)
    res <- merge(res, ann, by = "Gene_id", all = FALSE, sort = FALSE)

    # enforce column order: Gene_id, .Count, .Fpkm/.Cpm, annotation columns
    res <- res[c("Gene_id", paste0(sample_cols, ".Count"), norm_cols, ann_cols)]

    write.table(res, file = argv$output, sep = "\t", quote = FALSE, row.names = FALSE)
    message("[annotate] ", nrow(res), " genes x ", ncol(res), " cols -> ", argv$output)
}

main(argv)
