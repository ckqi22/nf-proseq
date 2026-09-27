#!/usr/bin/env Rscript
# =============================================================================
# pausing_index.R — Compute the PRO-seq pausing index for each transcript
#
# PI_i = (promoter_density + eps) / (genebody_density + eps)
#   promoter_density = promoter reads / promoter window length
#   genebody_density = genebody reads / genebody length
#   eps = 伪计数（分子分母同加，避免 genebody=0 时分母为 0）
#
# Input: two single-base count matrices produced by the pol2_profile subworkflow
#   (singlebase_count_merge.R), each with columns: transcript_id, gene_id, Length, <sample1>, ...
#   - promoter matrix:  'Length' = promoter window size (constant across genes)
#   - genebody matrix:  'Length' = genebody length (drives the min-length filter)
# Output: a single tab-separated PI table.
# =============================================================================
suppressWarnings(suppressMessages({
    library(argparser)
}))

# ── Parse args ──
argv <- arg_parser("Compute length-normalized PRO-seq pausing index")
argv <- add_argument(argv, "--promoter_count",  help = "Promoter count matrix (transcript_id, gene_id, Length, samples)")
argv <- add_argument(argv, "--genebody_count",  help = "Genebody count matrix (transcript_id, gene_id, Length, samples)")
argv <- add_argument(argv, "--min_genebody_len", help = "Minimum gene body length (bp)", default = 800, type = "integer")
argv <- add_argument(argv, "--pseudocount",     help = "Pseudocount added to both densities", default = 1e-3, type = "double")
argv <- add_argument(argv, "--output",          help = "Output PI table (tsv)")
argv <- parse_args(argv)


# ── Compute pausing index ──
compute_pausing_index <- function(promoter, genebody, min_len, eps) {
    gene_col <- colnames(promoter)[1]
    if (colnames(genebody)[1] != gene_col)
        stop("promoter and genebody count matrices must share the same first column (transcript_id)")
    if (!"Length" %in% colnames(promoter)) stop("promoter matrix must have a 'Length' column")
    if (!"Length" %in% colnames(genebody)) stop("genebody matrix must have a 'Length' column")

    # Align both matrices on common genes, then filter by gene body length.
    common <- intersect(promoter[[gene_col]], genebody[[gene_col]])
    promoter <- promoter[match(common, promoter[[gene_col]]), , drop = FALSE]
    genebody <- genebody[match(common, genebody[[gene_col]]), , drop = FALSE]

    genebody_len <- as.numeric(genebody[["Length"]])
    keep   <- genebody_len >= min_len
    promoter <- promoter[keep, , drop = FALSE]
    genebody <- genebody[keep, , drop = FALSE]

    genes   <- promoter[[gene_col]]
    promoter_len <- as.numeric(promoter[["Length"]])   # promoter window length
    genebody_len <- as.numeric(genebody[["Length"]])   # genebody length

    sample_cols <- setdiff(colnames(promoter), c(gene_col, "gene_id", "Length"))

    out <- data.frame(transcript_id = genes, gene_id = promoter[["gene_id"]],
                      promoter_len = promoter_len, genebody_len = genebody_len,
                      stringsAsFactors = FALSE)
    for (s in sample_cols) {
        pause <- as.numeric(promoter[[s]])
        body  <- as.numeric(genebody[[s]])
        # 分子分母同加伪计数 eps，避免 genebody=0 时分母为 0（PI 恒有定义）。
        pi <- (pause / promoter_len + eps) / (body / genebody_len + eps)
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

    out <- compute_pausing_index(promoter, genebody, min_len = as.integer(argv$min_genebody_len), eps = argv$pseudocount)

    write.table(out, file = argv$output, sep = "\t", quote = FALSE, row.names = FALSE)
    message("[pausing_index] ", nrow(out), " transcripts x ", ncol(out), " cols -> ", argv$output)
}

main(argv)
