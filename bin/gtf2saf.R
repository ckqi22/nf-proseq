#!/usr/bin/env Rscript
# =============================================================================
# gtf2saf.R — Generate per-gene genebody-union SAF from GTF
#
# 产物：
#   - genebody_union.saf : 每个 gene 所有 transcript 的 [TSS+offset, TES] 区间
#                          做 union（1-based 闭区间；同一 gene 可多行、GeneID 相同）。
#
# 供 featureCounts 定量。其余 BED 产物（tss/promoter/genebody/gene）
# 已移至 bin/gtf2bed.R。
#
# 坐标约定：
#   GTF / SAF = 1-based 闭区间。
#   正链 gene：TSS = start，TES = end；负链 gene：TSS = end，TES = start。
# =============================================================================
suppressWarnings(suppressMessages({
    library(argparser)
    library(rtracklayer)
    library(GenomicRanges)
}))

# ── Parse args ──
argv <- arg_parser("Generate genebody-union SAF from GTF")
argv <- add_argument(argv, "--gtf",             help = "Path to reference GTF file")
argv <- add_argument(argv, "--genebody_offset", help = "Gene body start offset (bp)", default = 301, type = "integer")
argv <- add_argument(argv, "--outdir",          help = "Output dir", default = "./")
argv <- parse_args(argv)

# ── Validate inputs ──
check_inputs <- function(argv) {
    if (!file.exists(argv$gtf)) stop("[gtf2saf] GTF file not found: ", argv$gtf)
}

# ── Read transcript records from the reference GTF ──
read_transcripts <- function(gtf_file) {
    message("[gtf2saf] Reading GTF: ", gtf_file)
    gr <- rtracklayer::import(gtf_file, format = "gtf")
    transcript_gr <- gr[gr$type == "transcript"]
    if (is.null(mcols(transcript_gr)$gene_id))
        stop("[gtf2saf] GTF missing 'gene_id' attribute on transcript entries")
    # 读取时即按坐标排序：sortSeqlevels 给自然序 seqlevels（chr1<...<chr10<chrX<chrM），
    # sort 再按 start 排 ranges（ignore.strand 保持纯坐标序，不按链分组）。
    transcript_gr <- sort(sortSeqlevels(transcript_gr), ignore.strand = TRUE)
    message("[gtf2saf] ", length(transcript_gr), " transcript entries")
    transcript_gr
}

# ── transcript 级清洗 ──
filter_transcripts <- function(gr, offset) {
    gene_id <- mcols(gr)$gene_id
    gr <- gr[!is.na(gene_id) & gene_id != "", ]      # 过滤空 gene_id
    gr[end(gr) - start(gr) >= offset, ]               # gene body 非空的最小长度（span >= offset）
}

# ── 所有 transcript 的 gene body 区间 union → SAF ──
build_genebody_union_saf <- function(gr, offset) {
    gene_id <- mcols(gr)$gene_id
    strand <- as.character(strand(gr))
    st <- start(gr); en <- end(gr)
    # 每个 transcript 的 gene body = [TSS+offset, TES]（1-based 闭区间）：
    #   + 链 [start+offset, end]；- 链 [start, end-offset]。排除 pause 区、含 intron。
    gb_start <- ifelse(strand == "+", st + offset, st)
    gb_end   <- ifelse(strand == "+", en,          en - offset)

    gb_gr <- GRanges(seqnames = seqnames(gr),
                     ranges   = IRanges(gb_start, gb_end),
                     strand   = strand)

    # 按 gene_id 分组后 reduce：同 gene 内合并重叠/相邻区间，跨 gene 不合并。
    gb_by_gene <- split(gb_gr, gene_id)
    union_list <- GenomicRanges::reduce(gb_by_gene)     # GRangesList，names = gene_id
    union_gr   <- unlist(union_list)

    # 每个 gene 一条链；从 transcript 记录取该 gene 的 strand，再回填到每个 reduced 区间。
    strand_by_gene <- strand[!duplicated(gene_id)]
    names(strand_by_gene) <- gene_id[!duplicated(gene_id)]

    n_per_gene  <- lengths(union_list)
    gene_id_rep <- rep(names(union_list), n_per_gene)
    strand_out  <- strand_by_gene[gene_id_rep]

    saf <- data.frame(
        GeneID = gene_id_rep,
        Chr    = as.character(seqnames(union_gr)),
        Start  = start(union_gr),                       # SAF：1-based 闭区间，直接写 GTF 坐标
        End    = end(union_gr),
        Strand = strand_out,
        stringsAsFactors = FALSE
    )
    saf <- saf[order(saf$Chr, saf$Start, saf$End), ]
    saf
}

# ── Main ──
main <- function(argv) {
    check_inputs(argv)

    transcript_gr <- filter_transcripts(read_transcripts(argv$gtf), argv$genebody_offset)

    union_saf <- build_genebody_union_saf(transcript_gr, argv$genebody_offset)
    write.table(union_saf, file = file.path(argv$outdir, "genebody_union.saf"),
                sep = "\t", quote = FALSE, row.names = FALSE)
    message("[gtf2saf] Genebody union SAF: ", nrow(union_saf), " intervals → ",
            file.path(argv$outdir, "genebody_union.saf"))

    message("[gtf2saf] Done.")
}

main(argv)
