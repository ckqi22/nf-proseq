#!/usr/bin/env Rscript
# =============================================================================
# featurecounts_merge.R — Merge per-sample featureCounts count files into one
# clean count matrix.
#
# Each input file (per-sample featureCounts output, SAF mode) has columns:
#   Geneid  Chr  Start  End  Strand  Length  <sample>.bam
# This merges them on Geneid into:
#   gene_id  length  <sample1>  <sample2> ...
#
#   - skip the leading '# Program:featureCounts ...' comment line(s)
#   - carry 'length' from the first file (annotation-derived, identical across
#     samples) so it stays a single column
#   - strip trailing '.bam' from sample columns so they match sample names
# =============================================================================
suppressWarnings(suppressMessages({
    library(argparser)
}))

# ── Parse args ──
argv <- arg_parser("Merge per-sample featureCounts files into a count matrix")
argv <- add_argument(argv, "--inputs", help = "Comma-separated per-sample featureCounts files")
argv <- add_argument(argv, "--output", help = "Output merged count matrix (tsv)")
argv <- parse_args(argv)

# ── Read one per-sample featureCounts file ──
read_counts <- function(f) {
    dt <- read.delim(f, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
                     comment.char = "#", check.names = FALSE)

    ann_cols <- c("Geneid", "Chr", "Start", "End", "Strand", "Length")
    cnt_cols <- setdiff(names(dt), ann_cols)
    if (length(cnt_cols) != 1)
        stop("[featurecounts_merge] expected 1 count column in ", f,
             ", got ", length(cnt_cols))

    # 先剥 R1-only BAM 的 ".r1.bam" 后缀，再剥普通 ".bam"，否则列名会残留
    # ".r1" 与 samplesheet 的样本名不匹配。
    sample <- sub("\\.r1\\.bam$", "", cnt_cols[1], ignore.case = TRUE)
    sample <- sub("\\.bam$", "", sample, ignore.case = TRUE)
    list(gene_id = dt$Geneid, length = dt$Length, sample = sample, count = dt[[cnt_cols[1]]])
}

# ── Main ──
main <- function(argv) {
    files <- trimws(unlist(strsplit(argv$inputs, ",", fixed = TRUE)))
    files <- files[nzchar(files)]
    if (length(files) == 0) stop("[featurecounts_merge] no input count files given")
    for (f in files) if (!file.exists(f)) stop("[featurecounts_merge] count file not found: ", f)

    counts_list <- lapply(files, read_counts)

    # Length is annotation-derived and identical across samples; take it once
    # from the first file and re-attach after merging.
    length_by_gene <- as.numeric(counts_list[[1]]$length)
    names(length_by_gene) <- counts_list[[1]]$gene_id

    # Full outer join on gene_id (defensive: files share the same SAF, so gene
    # sets are identical in practice, but all=TRUE keeps any stray gene).
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
    message("[featurecounts_merge] ", nrow(res), " genes x ", ncol(res),
            " cols -> ", argv$output)
}

main(argv)
