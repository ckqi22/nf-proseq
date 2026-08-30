#!/usr/bin/env Rscript
# =============================================================================
# singlebase_count_merge.R — Merge per-sample single-base count files into
# one count matrix (gene_id length samples). Generic over region type via --suffix.
#
# Each per-sample file (headerless, tab-separated) has columns:
#   gene_id  length  count
# Merges on gene_id into:
#   gene_id  length  <sample1>  <sample2> ...
#   - sample name is derived from the filename by stripping --suffix
#   - 'length' is carried from the first file (annotation-derived, identical
#     across samples since every sample uses the same region BED)
# =============================================================================
suppressWarnings(suppressMessages({
    library(argparser)
}))

# ── Parse args ──
argv <- arg_parser("Merge per-sample single-base region count files into a matrix")
argv <- add_argument(argv, "--inputs", help = "Comma-separated per-sample count files")
argv <- add_argument(argv, "--suffix", help = "Filename suffix to strip for sample name (e.g. .promoter.counts.txt)")
argv <- add_argument(argv, "--output", help = "Output merged count matrix (tsv)")
argv <- parse_args(argv)

# ── Read one per-sample count file ──
read_counts <- function(f, suffix) {
    dt <- read.delim(f, header = FALSE, sep = "\t", stringsAsFactors = FALSE,
                     col.names = c("gene_id", "length", "count"))
    suffix_esc <- gsub("\\.", "\\\\.", suffix)          # 把 . 转义为 \. 供正则匹配
    sample <- sub(paste0(suffix_esc, "$"), "", basename(f))
    list(gene_id = dt$gene_id, length = dt$length, sample = sample, count = dt$count)
}

# ── Main ──
main <- function(argv) {
    files <- trimws(unlist(strsplit(argv$inputs, ",", fixed = TRUE)))
    files <- files[nzchar(files)]
    if (length(files) == 0) stop("[singlebase_count_merge] no input count files given")
    for (f in files) if (!file.exists(f)) stop("[singlebase_count_merge] count file not found: ", f)

    counts_list <- lapply(files, read_counts, suffix = argv$suffix)

    # Length is annotation-derived and identical across samples; take it once
    # from the first file and re-attach after merging.
    length_by_gene <- as.numeric(counts_list[[1]]$length)
    names(length_by_gene) <- counts_list[[1]]$gene_id

    # Full outer join on gene_id (defensive: all samples share the same region
    # BED, so gene sets are identical in practice).
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
    message("[singlebase_count_merge] ", nrow(res), " genes x ", ncol(res),
            " cols -> ", argv$output)
}

main(argv)
