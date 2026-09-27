#!/usr/bin/env Rscript
# =============================================================================
# singlebase_count_merge.R — Merge per-sample single-base count files into
# one count matrix (transcript_id gene_id Length samples). Sample names are
# supplied explicitly via --sample (parallel to --count), not derived from
# filenames. gene_id is joined from --tx2gene (gtf2bed.R tx2gene.tsv).
#
# Each per-sample file (headerless, tab-separated) has columns:
#   transcript_id  length  count
# Merges on transcript_id into:
#   transcript_id  gene_id  Length  <sample1>  <sample2> ...
#   - sample names come from --sample, one per --count file, in order
#   - 'Length' is carried from the first file (annotation-derived, identical
#     across samples since every sample uses the same region BED)
#   - 'gene_id' is joined from --tx2gene (transcript_id → gene_id)
# =============================================================================
suppressWarnings(suppressMessages({
    library(argparser)
}))

# ── Parse args ──
argv <- arg_parser("Merge per-sample single-base region count files into a matrix")
argv <- add_argument(argv, "--sample",  help = "Comma-separated sample names (one per count file, in order)")
argv <- add_argument(argv, "--count",   help = "Comma-separated per-sample count files (same order as --sample)")
argv <- add_argument(argv, "--tx2gene", help = "transcript_id -> gene_id map (gtf2bed.R tx2gene.tsv, no header)")
argv <- add_argument(argv, "--output",  help = "Output merged count matrix (tsv)")
argv <- parse_args(argv)

# ── Read one per-sample count file ──
read_counts <- function(f) {
    dt <- read.delim(f, header = FALSE, sep = "\t", stringsAsFactors = FALSE,
                     col.names = c("transcript_id", "length", "count"))
    list(transcript_id = dt$transcript_id, length = dt$length, count = dt$count)
}

# ── Main ──
main <- function(argv) {
    samples <- trimws(unlist(strsplit(argv$sample, ",", fixed = TRUE)))
    samples <- samples[nzchar(samples)]
    files <- trimws(unlist(strsplit(argv$count, ",", fixed = TRUE)))
    files <- files[nzchar(files)]

    if (length(samples) == 0) stop("[singlebase_count_merge] no sample names given (--sample)")
    if (length(files) == 0)   stop("[singlebase_count_merge] no count files given (--count)")
    if (length(samples) != length(files))
        stop("[singlebase_count_merge] --sample (", length(samples), ") and --count (",
             length(files), ") must have the same length")
    for (f in files) if (!file.exists(f)) stop("[singlebase_count_merge] count file not found: ", f)

    counts_list <- lapply(files, read_counts)

    # Length is annotation-derived and identical across samples; take it once
    # from the first file and re-attach after merging.
    length_by_tx <- as.numeric(counts_list[[1]]$length)
    names(length_by_tx) <- counts_list[[1]]$transcript_id

    # Full outer join on transcript_id (defensive: all samples share the same region
    # BED, so transcript sets are identical in practice).
    res <- data.frame(transcript_id = counts_list[[1]]$transcript_id, stringsAsFactors = FALSE)
    for (i in 1:length(counts_list)) {
        p <- counts_list[[i]]
        one <- data.frame(transcript_id = p$transcript_id, count = p$count, stringsAsFactors = FALSE)
        names(one) <- c("transcript_id", samples[i])
        res <- merge(res, one, by = "transcript_id", all = TRUE, sort = FALSE)
    }

    res$Length <- length_by_tx[match(res$transcript_id, names(length_by_tx))]

    # gene_id 由 tx2gene 映射加入（transcript_id → gene_id）
    if (is.null(argv$tx2gene) || !file.exists(argv$tx2gene))
        stop("[singlebase_count_merge] --tx2gene required and must exist")
    txg <- read.delim(argv$tx2gene, header = FALSE, sep = "\t", stringsAsFactors = FALSE,
                      col.names = c("transcript_id", "gene_id"))
    txg <- txg[!duplicated(txg$transcript_id), , drop = FALSE]
    res$gene_id <- txg$gene_id[match(res$transcript_id, txg$transcript_id)]
    if (anyNA(res$gene_id))
        stop("[singlebase_count_merge] ", sum(is.na(res$gene_id)),
             " transcript_id(s) missing from tx2gene map")

    sample_cols <- setdiff(names(res), c("transcript_id", "gene_id", "Length"))
    res <- res[c("transcript_id", "gene_id", "Length", sample_cols)]

    write.table(res, file = argv$output, sep = "\t", quote = FALSE, row.names = FALSE)
    message("[singlebase_count_merge] ", nrow(res), " transcripts x ", ncol(res),
            " cols -> ", argv$output)
}

main(argv)
