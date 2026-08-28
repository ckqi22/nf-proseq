#!/usr/bin/env Rscript
# =============================================================================
# genebody_bedtools_merge.R — Merge per-sample bedtools single-base gene-body
# count files into one count matrix.
#
# Each per-sample file (headerless, tab-separated) has columns:
#   gene_id  length  count
# This merges them on gene_id into:
#   gene_id  length  <sample1>  <sample2> ...
#
#   - sample name is derived from the filename (strip ".genebody.counts.txt")
#   - 'length' is carried from the first file (annotation-derived, identical
#     across samples since every sample uses the same genebody.bed)
# =============================================================================
suppressWarnings(suppressMessages({
    library(argparser)
}))

# ── Parse args ──
argv <- arg_parser("Merge per-sample bedtools gene-body count files into a matrix")
argv <- add_argument(argv, "--inputs", help = "Comma-separated per-sample count files")
argv <- add_argument(argv, "--output", help = "Output merged count matrix (tsv)")
argv <- parse_args(argv)

# ── Read one per-sample count file ──
read_counts <- function(f) {
    dt <- read.delim(f, header = FALSE, sep = "\t", stringsAsFactors = FALSE,
                     col.names = c("gene_id", "length", "count"))
    sample <- sub("\\.genebody\\.counts\\.txt$", "", basename(f))
    list(gene_id = dt$gene_id, length = dt$length, sample = sample, count = dt$count)
}

# ── Main ──
main <- function(argv) {
    files <- trimws(unlist(strsplit(argv$inputs, ",", fixed = TRUE)))
    files <- files[nzchar(files)]
    if (length(files) == 0) stop("[genebody_bedtools_merge] no input count files given")
    for (f in files) if (!file.exists(f)) stop("[genebody_bedtools_merge] count file not found: ", f)

    counts_list <- lapply(files, read_counts)

    # Length is annotation-derived and identical across samples; take it once
    # from the first file and re-attach after merging.
    length_by_gene <- as.numeric(counts_list[[1]]$length)
    names(length_by_gene) <- counts_list[[1]]$gene_id

    # Full outer join on gene_id (defensive: all samples share the same
    # genebody.bed, so gene sets are identical in practice).
    res <- data.frame(gene_id = counts_list[[1]]$gene_id, stringsAsFactors = FALSE)
    for (p in counts_list) {
        one <- data.frame(gene_id = p$gene_id, count = p$count, stringsAsFactors = FALSE)
        names(one) <- c("gene_id", p$sample)
        res <- merge(res, one, by = "gene_id", all = TRUE, sort = FALSE)
    }

    res$length <- length_by_gene[match(res$gene_id, names(length_by_gene))]
    sample_cols <- setdiff(names(res), c("gene_id", "length"))
    res <- res[c("gene_id", "length", sample_cols)]

    write.table(res, file = argv$output, sep = "\t", quote = FALSE, row.names = FALSE)
    message("[genebody_bedtools_merge] ", nrow(res), " genes x ", ncol(res),
            " cols -> ", argv$output)
}

main(argv)
