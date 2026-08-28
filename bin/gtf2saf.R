#!/usr/bin/env Rscript
# =============================================================================
# gtf2saf.R — Generate promoter and genebody SAF annotations from GTF
#
# Reads GTF once, outputs two SAF files:
#   1. promoter window:      TSS-tss_upstream    to TSS+tss_downstream
#   2. genebody window:      TSS+genebody_offset to TES
#
# SAF Format (featureCounts requires these exact column names, case-sensitive)
# GeneID    Chr Start   End Strand
# GeneA chr1    100 400 +
# GeneB chr2    1500    2000    -
#
# Plus-strand genes:  TSS = start, TES = end
# Minus-strand genes: TSS = end,   TES = start
# =============================================================================
suppressWarnings(suppressMessages({
    library(argparser)
    library(rtracklayer)
}))

# ── Parse args ──
argv <- arg_parser("Generate promoter and genebody SAF annotations from GTF")
argv <- add_argument(argv, "--gtf",              help = "Path to reference GTF file")
argv <- add_argument(argv, "--tss_upstream",     help = "TSS upstream window (bp)",    default = 50,  type = "integer")
argv <- add_argument(argv, "--tss_downstream",   help = "TSS downstream window (bp)",  default = 300, type = "integer")
argv <- add_argument(argv, "--genebody_offset",  help = "Gene body start offset (bp)", default = 301, type = "integer")
argv <- add_argument(argv, "--outdir",           help = "Output dir", default = "./")
argv <- parse_args(argv)

# ── Read GTF ──
read_gtf <- function(gtf_file) {
    message("[gtf2saf] Reading GTF: ", gtf_file)
    gr <- rtracklayer::import(gtf_file, format = "gtf")   # GRanges; attributes parsed into mcols
    gr <- gr[gr$type == "gene"]                           # keep gene-level entries only
    if (is.null(mcols(gr)$gene_id)) stop("GTF missing 'gene_id' attribute")
    # 读取时即按坐标排序：sortSeqlevels 给自然序 seqlevels（chr1<...<chr10<chrX<chrY<chrM），
    # sort 再按 start 排 ranges（ignore.strand 保持纯坐标序，不按链分组）。
    gr <- sortSeqlevels(gr)
    gr <- sort(gr, ignore.strand = TRUE)
    message("[gtf2saf] ", length(gr), " gene entries")
    return(gr)
}

# ── Validate inputs ──
check_inputs <- function(argv) {
    if (!file.exists(argv$gtf)) stop("GTF file not found: ", argv$gtf)
    if (argv$genebody_offset <= argv$tss_downstream)
        stop("genebody_offset (", argv$genebody_offset, ") must be > tss_downstream (",
             argv$tss_downstream, ") so the genebody does not overlap the promoter window")
}

filter_gtf <- function(gtf, genebody_offset) {
    gene_len <- end(gtf) - start(gtf) + 1

    gene_id <- mcols(gtf)$gene_id
    ok_id   <- !is.na(gene_id) & gene_id != ""

    keep <- gene_len > genebody_offset + 1 & ok_id
    message("[gtf2saf] filter_gtf: removed ", sum(!keep),
            " genes (short gene_len <= ", genebody_offset + 1, " bp, or empty gene_id); ",
            sum(keep), " kept")
    gtf[keep, ]
}


# ── Build SAF data.frame from per-strand coordinate vectors ──
build_saf <- function(gene_ids, chr, strand, plus_start, plus_end, minus_start, minus_end) {
    plus <- strand == "+"
    data.frame(
        GeneID = gene_ids, Chr = chr,
        Start  = pmax(ifelse(plus, plus_start, minus_start), 1),
        End    = ifelse(plus, plus_end, minus_end),
        Strand = strand, stringsAsFactors = FALSE
    )
}

# ── Main ──
main <- function(argv) {
    check_inputs(argv)

    tss_upstream     <- argv$tss_upstream
    tss_downstream   <- argv$tss_downstream
    genebody_offset  <- argv$genebody_offset

    gr <- read_gtf(argv$gtf)
    gr <- filter_gtf(gr, genebody_offset)

    gene_ids <- mcols(gr)$gene_id
    chr    <- as.character(seqnames(gr))
    start  <- start(gr)
    end    <- end(gr)
    strand <- as.character(strand(gr))

    # ── Build TSS SAF ──
    tss_saf <- build_saf(
        gene_ids, chr, strand,
        plus_start  = start - tss_upstream,
        plus_end    = pmin(start + tss_downstream, end),   # clamp to gene end (TES)
        minus_start = pmax(end   - tss_downstream, start), # clamp to gene start (TES)
        minus_end   = end   + tss_upstream
    )

    # ── Build gene body SAF ──
    gb_saf <- build_saf(
        gene_ids, chr, strand,
        plus_start  = start + genebody_offset, plus_end  = end,
        minus_start = start,                    minus_end = end - genebody_offset
    )

    # ── Deduplicate ──
    tss_saf <- tss_saf[!duplicated(tss_saf$GeneID), ]
    gb_saf  <- gb_saf[!duplicated(gb_saf$GeneID), ]

    # ── Write outputs ──
    write.table(tss_saf, file = file.path(argv$outdir, "promoter.saf"), sep = "\t", quote = FALSE, row.names = FALSE)
    message("[gtf2saf] Promoter SAF: ", nrow(tss_saf), " regions → ", file.path(argv$outdir, "promoter.saf"))

    write.table(gb_saf, file = file.path(argv$outdir, "genebody.saf"), sep = "\t", quote = FALSE, row.names = FALSE)
    message("[gtf2saf] Genebody SAF: ", nrow(gb_saf), " regions → ", file.path(argv$outdir, "genebody.saf"))

    # ── Gene body BED（0-based 半开，供 bedtools map 与 genomecov -bg 对齐）──
    #   SAF 为 1-based 闭区间，BED 为 0-based 半开：Start = saf_start - 1, End = saf_end。
    #   基因体长度保持 End - Start + 1，与 featureCounts 的 Length 一致。
    gb_bed <- data.frame(
        Chr = gb_saf$Chr, Start = gb_saf$Start - 1, End = gb_saf$End,
        GeneID = gb_saf$GeneID, Score = ".", Strand = gb_saf$Strand,
        stringsAsFactors = FALSE
    )
    write.table(gb_bed, file = file.path(argv$outdir, "genebody.bed"), sep = "\t",
                quote = FALSE, row.names = FALSE, col.names = FALSE)
    message("[gtf2saf] Genebody BED: ", nrow(gb_bed), " regions → ", file.path(argv$outdir, "genebody.bed"))

    # ── Full-gene BED（0-based 半开，供 POL2_COUNT 单碱基 5' 计数）──
    #   全基因跨度（start→end），BED6 列序：Chr Start End GeneID . Strand。
    #   基因集与 genebody.bed 一致（同为 filter_gtf 之后），保证下游 gene_id 全集一致。
    gene_bed <- data.frame(
        Chr = chr, Start = start - 1, End = end,
        GeneID = gene_ids, Score = ".", Strand = strand,
        stringsAsFactors = FALSE
    )
    write.table(gene_bed, file = file.path(argv$outdir, "gene.bed"), sep = "\t",
                quote = FALSE, row.names = FALSE, col.names = FALSE)
    message("[gtf2saf] Gene BED: ", nrow(gene_bed), " regions → ", file.path(argv$outdir, "gene.bed"))

    # ── Strand-split gene BED（供 tss_meta 链特异性 computeMatrix）──
    #   PRO-seq reverse-stranded: + 基因信号在 plus.bigWig, - 基因信号在 minus.bigWig。
    #   tss_meta 需要按链拆分基因，分别用正确的 bigWig 做 metagene。
    plus_idx  <- gene_bed$Strand == "+"
    minus_idx <- gene_bed$Strand == "-"
    write.table(gene_bed[plus_idx, ],  file = file.path(argv$outdir, "plus_genes.bed"), sep = "\t",
                quote = FALSE, row.names = FALSE, col.names = FALSE)
    write.table(gene_bed[minus_idx, ], file = file.path(argv$outdir, "minus_genes.bed"), sep = "\t",
                quote = FALSE, row.names = FALSE, col.names = FALSE)
    message("[gtf2saf] Strand-split BED: + ", sum(plus_idx), " genes, - ", sum(minus_idx), " genes")

    message("[gtf2saf] Done.")
}

main(argv)
