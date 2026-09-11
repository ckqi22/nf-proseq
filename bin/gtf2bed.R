#!/usr/bin/env Rscript
# =============================================================================
# gtf2bed.R — Generate TSS / promoter / genebody / gene BED from GTF
#
# 产出（均为 BED6，0-based 半开）：
#   - tss.bed       : 每个代表 transcript 的 TSS 碱基（name = gene_id）
#   - promoter.bed  : TSS - tss_upstream  → TSS + tss_downstream（代表 transcript）
#   - genebody.bed  : TSS + genebody_offset → TES（代表 transcript）
#   - gene.bed      : 每个 gene 的 [start, end] 跨度（参考 GTF 的 gene 记录）
#
# 其中 TSS/TES 取 --rep_gtf（longest_tx.R 产出）的代表 transcript 坐标。
#
# 坐标约定：
#   GTF = 1-based 闭区间；BED = 0-based 半开（Start = 1-based_start - 1）。
#   正链 gene：TSS = start，TES = end；负链 gene：TSS = end，TES = start。
# =============================================================================
suppressWarnings(suppressMessages({
    library(argparser)
    library(rtracklayer)
    library(GenomicRanges)
}))

# ── Parse args ──
argv <- arg_parser("Generate TSS/promoter/genebody/gene BED from GTF")
argv <- add_argument(argv, "--gtf",             help = "Path to reference GTF file")
argv <- add_argument(argv, "--rep_gtf",         help = "Representative transcript GTF (from longest_tx.R)")
argv <- add_argument(argv, "--tss_upstream",    help = "TSS upstream window (bp)",    default = 50,  type = "integer")
argv <- add_argument(argv, "--tss_downstream",  help = "TSS downstream window (bp)",  default = 300, type = "integer")
argv <- add_argument(argv, "--genebody_offset", help = "Gene body start offset (bp)", default = 301, type = "integer")
argv <- add_argument(argv, "--outdir",          help = "Output dir", default = "./")
argv <- parse_args(argv)

# ── Validate inputs ──
check_inputs <- function(argv) {
    if (!file.exists(argv$gtf)) stop("[gtf2bed] GTF file not found: ", argv$gtf)
    if (is.null(argv$rep_gtf) || !file.exists(argv$rep_gtf))
        stop("[gtf2bed] Representative transcript GTF not found: ", argv$rep_gtf)
    if (argv$genebody_offset <= argv$tss_downstream)
        stop("[gtf2bed] genebody_offset (", argv$genebody_offset, ") must be > tss_downstream (",
             argv$tss_downstream, ") so the genebody does not overlap the promoter window")
}

# ── Read gene records from the reference GTF ──
read_genes <- function(gtf_file) {
    message("[gtf2bed] Reading reference GTF (gene records): ", gtf_file)
    gr <- rtracklayer::import(gtf_file, format = "gtf")
    gene_gr <- gr[gr$type == "gene"]
    if (is.null(mcols(gene_gr)$gene_id))
        stop("[gtf2bed] GTF missing 'gene_id' attribute on gene entries")
    gene_gr <- sort(sortSeqlevels(gene_gr), ignore.strand = TRUE)
    message("[gtf2bed] ", length(gene_gr), " gene entries")
    gene_gr
}

# ── Read representative transcript GTF (longest_tx.R 产物)，返回全部 transcript ──
read_rep_tx <- function(rep_gtf) {
    message("[gtf2bed] Reading representative transcript GTF: ", rep_gtf)
    gr <- rtracklayer::import(rep_gtf, format = "gtf")
    tx <- gr[gr$type == "transcript"]
    if (length(tx) == 0)
        stop("[gtf2bed] no 'transcript' entries in representative GTF: ", rep_gtf)
    if (is.null(mcols(tx)$gene_id))
        stop("[gtf2bed] representative GTF missing 'gene_id' on transcript entries")
    if (is.null(mcols(tx)$transcript_id))
        stop("[gtf2bed] representative GTF missing 'transcript_id' on transcript entries")
    message("[gtf2bed] ", length(tx), " representative transcripts")
    tx
}

# ── Build BED data.frame from per-strand coordinate vectors（1-based 闭区间）──
build_bed <- function(gene_ids, chr, strand, plus_start, plus_end, minus_start, minus_end) {
    plus <- strand == "+"
    data.frame(
        Chr    = chr,
        Start  = pmax(ifelse(plus, plus_start, minus_start), 1) - 1,   # 1-based 闭 → 0-based 半开
        End    = ifelse(plus, plus_end, minus_end),
        GeneID = gene_ids,
        Score  = 0,
        Strand = strand,
        stringsAsFactors = FALSE
    )
}

# ── Main ──
main <- function(argv) {
    check_inputs(argv)

    gene_gr <- read_genes(argv$gtf)

    # ── 代表 transcript（全部）→ tss.bed ──
    rep_tx  <- read_rep_tx(argv$rep_gtf)
    gene_ids <- mcols(rep_tx)$gene_id
    chr      <- as.character(seqnames(rep_tx))
    strand   <- as.character(strand(rep_tx))
    start    <- start(rep_tx)
    end      <- end(rep_tx)

    # TSS 碱基（1-based）：+ = start；- = end；BED 半开 [tss-1, tss)
    tss <- ifelse(strand == "-", end, start)
    tss_bed <- data.frame(
        Chr = chr, Start = tss - 1, End = tss,
        GeneID = gene_ids, Score = 0, Strand = strand,
        stringsAsFactors = FALSE
    )
    write.table(tss_bed, file = file.path(argv$outdir, "tss.bed"), sep = "\t",
                quote = FALSE, row.names = FALSE, col.names = FALSE)
    message("[gtf2bed] TSS BED: ", nrow(tss_bed), " regions → ",
            file.path(argv$outdir, "tss.bed"))

    # ── 代表 transcript（offset 过滤后）→ promoter.bed / genebody.bed ──
    keep    <- (end - start) >= argv$genebody_offset   # gene body 非空（span >= offset）
    gids    <- gene_ids[keep]
    gchr    <- chr[keep]
    gstrand <- strand[keep]
    gstart  <- start[keep]
    gend    <- end[keep]

    # promoter 窗口（1-based 闭区间）：+ [start-50, start+300]；- [end-300, end+50]，clamp 到 transcript 边界。
    prom_bed <- build_bed(
        gids, gchr, gstrand,
        plus_start  = gstart - argv$tss_upstream,
        plus_end    = pmin(gstart + argv$tss_downstream, gend),
        minus_start = pmax(gend   - argv$tss_downstream, gstart),
        minus_end   = gend   + argv$tss_upstream
    )
    # gene body 窗口（1-based 闭区间）：+ [start+301, end]；- [start, end-301]。
    gb_bed <- build_bed(
        gids, gchr, gstrand,
        plus_start  = gstart + argv$genebody_offset, plus_end  = gend,
        minus_start = gstart,                        minus_end = gend - argv$genebody_offset
    )

    write.table(prom_bed, file = file.path(argv$outdir, "promoter.bed"), sep = "\t",
                quote = FALSE, row.names = FALSE, col.names = FALSE)
    message("[gtf2bed] Promoter BED: ", nrow(prom_bed), " regions → ",
            file.path(argv$outdir, "promoter.bed"))

    write.table(gb_bed, file = file.path(argv$outdir, "genebody.bed"), sep = "\t",
                quote = FALSE, row.names = FALSE, col.names = FALSE)
    message("[gtf2bed] Genebody BED: ", nrow(gb_bed), " regions → ",
            file.path(argv$outdir, "genebody.bed"))

    # ── gene 级产物（gene.bed，供信号表 intersect）──
    gene_bed <- data.frame(
        Chr = as.character(seqnames(gene_gr)), Start = start(gene_gr) - 1, End = end(gene_gr),
        GeneID = mcols(gene_gr)$gene_id, Score = 0, Strand = as.character(strand(gene_gr)),
        stringsAsFactors = FALSE
    )
    write.table(gene_bed, file = file.path(argv$outdir, "gene.bed"), sep = "\t",
                quote = FALSE, row.names = FALSE, col.names = FALSE)
    message("[gtf2bed] Gene BED: ", nrow(gene_bed), " regions → ",
            file.path(argv$outdir, "gene.bed"))

    message("[gtf2bed] Done.")
}

main(argv)
