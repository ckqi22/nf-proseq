#!/usr/bin/env Rscript
# =============================================================================
# longest_tx.R — Select one representative transcript
# per gene from a GENCODE GTF (longest transcript by summed exon length),
# and write a GTF containing only those representative transcripts.
#
# 每个 gene 选一个「代表 transcript」作为 TSS 参考：
#   - 从 GTF 筛选指定 gene_type（默认 protein_coding）；
#   - 每个 transcript 的长度 = 其所有 exon 长度之和；
#   - 每个 gene 选 exon 总长最长的 transcript；
#   - 多个 transcript 长度相同时保留稳定排序后的第一个（order(gene_id, -tx_len) + !duplicated，即 GTF 文件顺序）；
#   - + 链 TSS = transcript start；- 链 TSS = transcript end。
#
# Input:
#   --gtf        参考 GTF（可 gz）。必填。
#   --gene-type  gene biotype 过滤（GTF 属性 gene_type，回退 gene_biotype）。默认 protein_coding。
#   --outdir     输出目录。默认 ./。
#
# Output:
#   longest_tx.gtf  仅含代表 transcript 的 GTF
#      （每个 gene 保留 gene + 代表 transcript + 其 exon 记录）
#
# Usage:
#   Rscript longest_tx.R --gtf gencode.v46.annotation.gtf.gz --outdir ./
# =============================================================================
suppressWarnings(suppressMessages({
    library(argparser)
    library(rtracklayer)
    library(GenomicRanges)
}))

# ── Parse args ──
argv <- arg_parser("Select one representative transcript per gene from a GENCODE GTF")
argv <- add_argument(argv, "--gtf",       help = "Reference GTF (gzipped ok). Required.")
argv <- add_argument(argv, "--gene-type", help = "Gene biotype filter (gene_type, fallback gene_biotype)", default = "protein_coding")
argv <- add_argument(argv, "--outdir",    help = "Output directory", default = "./")
argv <- parse_args(argv)

# ── Main ──
main <- function(argv) {
    if (is.null(argv$gtf)) stop("[longest_tx] --gtf is required")
    if (!file.exists(argv$gtf)) stop("[longest_tx] GTF not found: ", argv$gtf)

    message("[longest_tx] Reading GTF: ", argv$gtf)
    gtf_gr <- rtracklayer::import(argv$gtf, format = "gtf")   # 按 GTF 文件顺序，不排序（平局依赖此顺序）

    # gene 记录 → gene_type map（GENCODE 用 gene_type；老版本用 gene_biotype）
    gene_rec <- gtf_gr[gtf_gr$type == "gene"]
    if (is.null(mcols(gene_rec)$gene_id))
        stop("[longest_tx] GTF missing 'gene_id' attribute on gene entries")
    gtype <- mcols(gene_rec)$gene_type
    if (is.null(gtype)) gtype <- mcols(gene_rec)$gene_biotype
    if (is.null(gtype))
        stop("[longest_tx] GTF missing 'gene_type'/'gene_biotype' on gene entries")
    gene_type_map <- setNames(as.character(gtype), mcols(gene_rec)$gene_id)

    # transcript 记录按 biotype 过滤（保持文件顺序）
    tx_rec <- gtf_gr[gtf_gr$type == "transcript"]
    if (is.null(mcols(tx_rec)$transcript_id))
        stop("[longest_tx] GTF missing 'transcript_id' attribute on transcript entries")
    tx_gtype <- gene_type_map[mcols(tx_rec)$gene_id]
    tx_rec <- tx_rec[!is.na(tx_gtype) & tx_gtype == argv$gene_type]
    if (length(tx_rec) == 0)
        stop("[longest_tx] no transcripts of gene type ", argv$gene_type, " found")

    # exon 长度 → 每个 transcript 的 exon 总长 → 每个 gene 选最长
    ex_rec <- gtf_gr[gtf_gr$type == "exon"]
    ex_len <- tapply(width(ex_rec), mcols(ex_rec)$transcript_id, sum)
    tx_len <- as.numeric(ex_len[mcols(tx_rec)$transcript_id])

    tx_df <- data.frame(gene_id           = as.character(mcols(tx_rec)$gene_id),
                        transcript_id     = as.character(mcols(tx_rec)$transcript_id),
                        transcript_length = tx_len,
                        stringsAsFactors  = FALSE)
    # 最长 transcript per gene（平局取稳定排序后第一个，即 GTF 文件顺序，不新增 tie-break key）
    tx_df <- tx_df[order(tx_df$gene_id, -tx_df$transcript_length), ]
    rep_df <- tx_df[!duplicated(tx_df$gene_id), ]

    # subset GTF：每个 gene 保留 gene + 代表 transcript + 其 exon 记录（保持原文件顺序）
    rep_gene_ids <- rep_df$gene_id
    rep_tx_ids   <- rep_df$transcript_id
    keep <- (gtf_gr$type == "gene"       & gtf_gr$gene_id       %in% rep_gene_ids) |
            (gtf_gr$type == "transcript" & gtf_gr$transcript_id %in% rep_tx_ids)   |
            (gtf_gr$type == "exon"       & gtf_gr$transcript_id %in% rep_tx_ids)
    out_gr <- gtf_gr[keep]

    out_file <- file.path(argv$outdir, "longest_tx.gtf")
    rtracklayer::export(out_gr, out_file, format = "gtf")
    message("[longest_tx] ", nrow(rep_df), " genes -> ", out_file)
    message("[longest_tx] Done.")
}

main(argv)
