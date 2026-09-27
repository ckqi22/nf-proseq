#!/usr/bin/env Rscript
# =============================================================================
# gtf2bed.R — Generate TSS / promoter / genebody BED + tx2gene map from GTF
#
# 产出（BED6 为 0-based 半开，name = transcript_id；tx2gene 为两列 TSV）：
#   - tss.bed      : 每个 protein_coding transcript 的 TSS 碱基
#   - promoter.bed : TSS - tss_upstream → TSS + tss_downstream
#   - genebody.bed : TSS + genebody_offset → TES
#   - tx2gene.tsv  : transcript_id → gene_id 映射（无表头两列）
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
argv <- arg_parser("Generate TSS/promoter/genebody BED + tx2gene map from GTF")
argv <- add_argument(argv, "--gtf",             help = "Path to reference GTF file")
argv <- add_argument(argv, "--tss_upstream",    help = "TSS upstream window (bp)",    default = 50,  type = "integer")
argv <- add_argument(argv, "--tss_downstream",  help = "TSS downstream window (bp)",  default = 300, type = "integer")
argv <- add_argument(argv, "--genebody_offset", help = "Gene body start offset (bp)", default = 301, type = "integer")
argv <- add_argument(argv, "--outdir",          help = "Output dir", default = "./")
argv <- parse_args(argv)

# ── Validate inputs ──
check_inputs <- function(argv) {
    if (!file.exists(argv$gtf)) stop("[gtf2bed] GTF file not found: ", argv$gtf)
    if (argv$genebody_offset <= argv$tss_downstream)
        stop("[gtf2bed] genebody_offset (", argv$genebody_offset, ") must be > tss_downstream (",
             argv$tss_downstream, ") so the genebody does not overlap the promoter window")
}

# ── Read all protein_coding transcripts ──
read_tx <- function(gtf_file) {
    message("[gtf2bed] Reading GTF: ", gtf_file)
    gr <- rtracklayer::import(gtf_file, format = "gtf")

    # gene 记录 → gene_type map（GENCODE 用 gene_type；老版本回退 gene_biotype）
    gene_rec <- gr[gr$type == "gene"]
    if (is.null(mcols(gene_rec)$gene_id))
        stop("[gtf2bed] GTF missing 'gene_id' attribute on gene entries")
    gtype <- mcols(gene_rec)$gene_type
    if (is.null(gtype)) gtype <- mcols(gene_rec)$gene_biotype
    if (is.null(gtype))
        stop("[gtf2bed] GTF missing 'gene_type'/'gene_biotype' on gene entries")
    gene_type_map <- setNames(as.character(gtype), mcols(gene_rec)$gene_id)

    # transcript 记录按 biotype 过滤
    tx <- gr[gr$type == "transcript"]
    if (length(tx) == 0)
        stop("[gtf2bed] no 'transcript' entries in GTF: ", gtf_file)
    if (is.null(mcols(tx)$transcript_id))
        stop("[gtf2bed] GTF missing 'transcript_id' attribute on transcript entries")
    if (is.null(mcols(tx)$gene_id))
        stop("[gtf2bed] GTF missing 'gene_id' attribute on transcript entries")
    tx_gtype <- gene_type_map[mcols(tx)$gene_id]
    tx <- tx[!is.na(tx_gtype) & tx_gtype == "protein_coding"]
    if (length(tx) == 0)
        stop("[gtf2bed] no protein_coding transcripts found")
    message("[gtf2bed] ", length(tx), " protein_coding transcripts")
    tx
}

# ── Build BED data.frame from per-strand coordinate vectors（1-based 闭区间）──
build_bed <- function(tx_ids, chr, strand, plus_start, plus_end, minus_start, minus_end) {
    plus <- strand == "+"
    data.frame(
        Chr    = chr,
        Start  = pmax(ifelse(plus, plus_start, minus_start), 1) - 1,   # 1-based 闭 → 0-based 半开
        End    = ifelse(plus, plus_end, minus_end),
        GeneID = tx_ids,
        Score  = 0,
        Strand = strand,
        stringsAsFactors = FALSE
    )
}

# ── Main ──
main <- function(argv) {
    check_inputs(argv)

    tx <- read_tx(argv$gtf)

    tx_ids  <- as.character(mcols(tx)$transcript_id)
    gene_ids <- as.character(mcols(tx)$gene_id)
    chr     <- as.character(seqnames(tx))
    strand  <- as.character(strand(tx))
    start   <- start(tx)
    end     <- end(tx)

    # TSS 碱基（1-based）：+ = start；- = end；BED 半开 [tss-1, tss)
    tss <- ifelse(strand == "-", end, start)
    tss_bed <- data.frame(
        Chr = chr, Start = tss - 1, End = tss,
        GeneID = tx_ids, Score = 0, Strand = strand,
        stringsAsFactors = FALSE
    )
    write.table(tss_bed, file = file.path(argv$outdir, "tss.bed"), sep = "\t",
                quote = FALSE, row.names = FALSE, col.names = FALSE)
    message("[gtf2bed] TSS BED: ", nrow(tss_bed), " regions → ",
            file.path(argv$outdir, "tss.bed"))

    # ── promoter / genebody（span >= offset 过滤）──
    keep    <- (end - start) >= argv$genebody_offset   # gene body 非空（span >= offset）
    txids   <- tx_ids[keep]
    gchr    <- chr[keep]
    gstrand <- strand[keep]
    gstart  <- start[keep]
    gend    <- end[keep]

    # promoter 窗口（1-based 闭区间）：+ [start-50, start+300]；- [end-300, end+50]，clamp 到 transcript 边界。
    prom_bed <- build_bed(
        txids, gchr, gstrand,
        plus_start  = gstart - argv$tss_upstream,
        plus_end    = pmin(gstart + argv$tss_downstream, gend),
        minus_start = pmax(gend   - argv$tss_downstream, gstart),
        minus_end   = gend   + argv$tss_upstream
    )
    # gene body 窗口（1-based 闭区间）：+ [start+301, end]；- [start, end-301]。
    gb_bed <- build_bed(
        txids, gchr, gstrand,
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

    # ── tx2gene.tsv：transcript_id → gene_id（无表头两列，按 transcript_id 去重）──
    tx2gene <- data.frame(transcript_id = tx_ids, gene_id = gene_ids, stringsAsFactors = FALSE)
    tx2gene <- tx2gene[!duplicated(tx2gene$transcript_id), , drop = FALSE]
    write.table(tx2gene, file = file.path(argv$outdir, "tx2gene.tsv"), sep = "\t",
                quote = FALSE, row.names = FALSE, col.names = FALSE)
    message("[gtf2bed] tx2gene map: ", nrow(tx2gene), " transcripts → ",
            file.path(argv$outdir, "tx2gene.tsv"))

    message("[gtf2bed] Done.")
}

main(argv)
