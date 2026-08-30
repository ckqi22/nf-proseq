#!/usr/bin/env Rscript
# =============================================================================
# gtf2saf.R — Generate per-gene promoter/genebody BED + genebody-union SAF from GTF
#
# 两分支产物（互不混淆）：
#   分支 A（单碱基定量 / PI，最长 transcript）：
#     - promoter.bed   : TSS - tss_upstream  → TSS + tss_downstream（0-based 半开）
#     - genebody.bed   : TSS + genebody_offset → TES（0-based 半开）
#     其中 TSS/TES 取该 gene 最长 transcript 的坐标。
#   分支 B（featureCounts 定量，所有 transcript 的 gene body union）：
#     - genebody_union.saf : 每个 gene 所有 transcript 的 [TSS+offset, TES] 区间
#                            做 union（1-based 闭区间；同一 gene 可多行、GeneID 相同）。
#
# 另有 gene 级产物（gene.bed，供 tss_meta 与逐碱基信号表 SIGNAL_TABLE intersect）：
#     - gene.bed
#
# 坐标约定：
#   GTF / SAF = 1-based 闭区间；BED = 0-based 半开（Start = 1-based_start - 1）。
#   正链 gene：TSS = start，TES = end；负链 gene：TSS = end，TES = start。
# =============================================================================
suppressWarnings(suppressMessages({
    library(argparser)
    library(rtracklayer)
    library(GenomicRanges)
}))

# ── Parse args ──
argv <- arg_parser("Generate promoter/genebody BED and genebody-union SAF from GTF")
argv <- add_argument(argv, "--gtf",             help = "Path to reference GTF file")
argv <- add_argument(argv, "--tss_upstream",    help = "TSS upstream window (bp)",    default = 50,  type = "integer")
argv <- add_argument(argv, "--tss_downstream",  help = "TSS downstream window (bp)",  default = 300, type = "integer")
argv <- add_argument(argv, "--genebody_offset", help = "Gene body start offset (bp)", default = 301, type = "integer")
argv <- add_argument(argv, "--outdir",          help = "Output dir", default = "./")
argv <- parse_args(argv)

# —— test args ——
setwd("/home/ck/bioinfo/proseq_test")
argv$gtf <- "/home/ck/bioinfo/proseq_test/Homo_sapiens.GRCh38.92_UCSC.gtf"

# ── Validate inputs ──
check_inputs <- function(argv) {
    if (!file.exists(argv$gtf)) stop("GTF file not found: ", argv$gtf)
    if (argv$genebody_offset <= argv$tss_downstream)
        stop("genebody_offset (", argv$genebody_offset, ") must be > tss_downstream (",
             argv$tss_downstream, ") so the genebody does not overlap the promoter window")
}

# ── Read GTF：拆出 gene 与 transcript 两类记录 ──
read_gtf <- function(gtf_file) {
    message("[gtf2saf] Reading GTF: ", gtf_file)
    gr <- rtracklayer::import(gtf_file, format = "gtf")   # GRanges；attributes 解析进 mcols
    gene_gr       <- gr[gr$type == "gene"]
    transcript_gr <- gr[gr$type == "transcript"]
    if (is.null(mcols(gene_gr)$gene_id))       stop("GTF missing 'gene_id' attribute on gene entries")
    if (is.null(mcols(transcript_gr)$gene_id)) stop("GTF missing 'gene_id' attribute on transcript entries")
    # 读取时即按坐标排序：sortSeqlevels 给自然序 seqlevels（chr1<...<chr10<chrX<chrM），
    # sort 再按 start 排 ranges（ignore.strand 保持纯坐标序，不按链分组）。
    gene_gr       <- sort(sortSeqlevels(gene_gr),       ignore.strand = TRUE)
    transcript_gr <- sort(sortSeqlevels(transcript_gr), ignore.strand = TRUE)
    message("[gtf2saf] ", length(gene_gr), " gene entries, ",
            length(transcript_gr), " transcript entries")
    list(gene = gene_gr, transcript = transcript_gr)
}

# ── transcript 级清洗 + 长度预计算（分支 A/B 共用，一次完成）──
filter_transcripts <- function(gr, offset) {
    gene_id <- mcols(gr)$gene_id
    gr <- gr[!is.na(gene_id) & gene_id != "", ]      # 过滤空 gene_id
    span <- end(gr) - start(gr)
    mcols(gr)$span <- span                            # 预计算长度，供下游复用
    gr[span >= offset, ]                              # gene body 非空的最小长度（span >= offset）
}

# ── 分支 A：每个 gene 选最长 transcript（相同长度则保留更上游的）──
select_longest_transcript <- function(gr) {
    gene_id <- mcols(gr)$gene_id
    tid     <- mcols(gr)$transcript_id
    strand  <- as.character(strand(gr))
    start   <- start(gr)
    end     <- end(gr)
    span    <- mcols(gr)$span       # 长度已在 filter_transcripts 里预计算
    # tie-break key（最上游 TSS）：+ 链最上游 = 最小 start，- 链最上游 = 最大 end。
    # 用 ifelse 造一个"同向升序"key，使统一升序即代表"最上游优先"。
    tss_key <- ifelse(strand == "+", start, -end)
    # 缺 transcript_id 时按坐标合成稳定 id，保证排序确定。
    if (is.null(tid) || all(is.na(tid))) {
        tid <- paste0(as.character(seqnames(gr)), ":", start, "-", end)
    } else {
        tid[is.na(tid)] <- paste0(as.character(seqnames(gr))[is.na(tid)], ":",
                                  start[is.na(tid)], "-", end[is.na(tid)])
    }

    df <- data.frame(i = seq_along(gr), gene_id = gene_id, span = span,
                     tss_key = tss_key, tid = tid, stringsAsFactors = FALSE)
    # 组内排序：span 降序 → tss_key 升序（最上游 TSS）→ transcript_id 字典序升序，取第一条。
    df <- df[order(df$gene_id, -df$span, df$tss_key, df$tid), ]
    keep <- !duplicated(df$gene_id)
    out <- gr[df$i[keep], ]
    mcols(out)$transcript_id <- df$tid[keep]   # 回填最终 transcript_id（含合成的），供检查输出
    out
}

# ── 分支 B：所有 transcript 的 gene body 区间 union → SAF ──
build_genebody_union_saf <- function(gr, offset) {
    gene_id <- mcols(gr)$gene_id
    strand <- as.character(strand(gr))
    st <- start(gr); en <- end(gr)
    # 每个 transcript 的 gene body = [TSS+offset, TES]（1-based 闭区间）：
    #   + 链 [start+offset, end]；- 链 [start, end-offset]。排除 pause 区、含 intron。
    # （空 gene_id 与短 transcript 已在 filter_transcripts 里过滤，此处无需再过滤。）
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

# ── Build SAF data.frame from per-strand coordinate vectors（1-based 闭区间）──
build_saf <- function(gene_ids, chr, strand, plus_start, plus_end, minus_start, minus_end) {
    plus <- strand == "+"
    data.frame(
        GeneID = gene_ids, Chr = chr,
        Start  = pmax(ifelse(plus, plus_start, minus_start), 1),
        End    = ifelse(plus, plus_end, minus_end),
        Strand = strand, stringsAsFactors = FALSE
    )
}

# ── Convert SAF (1-based closed) to BED (0-based half-open) ──
to_bed <- function(saf) data.frame(
    Chr = saf$Chr, Start = saf$Start - 1, End = saf$End,
    GeneID = saf$GeneID, Score = ".", Strand = saf$Strand,
    stringsAsFactors = FALSE
)

# ── Main ──
main <- function(argv) {
    check_inputs(argv)

    tss_upstream    <- argv$tss_upstream
    tss_downstream  <- argv$tss_downstream
    genebody_offset <- argv$genebody_offset

    g <- read_gtf(argv$gtf)
    gene_gr       <- g$gene                                              # gene 级：不再按长度过滤（gene.bed 保留所有 gene）
    transcript_gr <- filter_transcripts(g$transcript, genebody_offset)   # 清洗 + 长度预计算 + 短 transcript 过滤

    # ── 分支 A：最长 transcript → TSS → promoter.bed / genebody.bed ──
    rep_tr <- select_longest_transcript(transcript_gr)

    # 输出挑选的代表（最长）transcript，供人工检查（span = end - start，排序长度量度）
    rep_tab <- data.frame(
        gene_id       = mcols(rep_tr)$gene_id,
        transcript_id = mcols(rep_tr)$transcript_id,
        chr           = as.character(seqnames(rep_tr)),
        start         = start(rep_tr),
        end           = end(rep_tr),
        strand        = as.character(strand(rep_tr)),
        span          = mcols(rep_tr)$span,
        stringsAsFactors = FALSE
    )
    write.table(rep_tab, file = file.path(argv$outdir, "rep_transcript.tsv"), sep = "\t",
                quote = FALSE, row.names = FALSE)
    message("[gtf2saf] Representative transcript: ", nrow(rep_tab), " genes → ",
            file.path(argv$outdir, "rep_transcript.tsv"))

    gene_ids <- mcols(rep_tr)$gene_id
    chr      <- as.character(seqnames(rep_tr))
    strand   <- as.character(strand(rep_tr))
    start    <- start(rep_tr)
    end      <- end(rep_tr)

    # promoter 窗口（1-based 闭区间）：+ [start-50, start+300]；- [end-300, end+50]，clamp 到 transcript 边界。
    prom_saf <- build_saf(
        gene_ids, chr, strand,
        plus_start  = start - tss_upstream,
        plus_end    = pmin(start + tss_downstream, end),
        minus_start = pmax(end   - tss_downstream, start),
        minus_end   = end   + tss_upstream
    )
    # gene body 窗口（1-based 闭区间）：+ [start+301, end]；- [start, end-301]。
    gb_saf <- build_saf(
        gene_ids, chr, strand,
        plus_start  = start + genebody_offset, plus_end  = end,
        minus_start = start,                   minus_end = end - genebody_offset
    )

    # 转 0-based 半开 BED（Start = saf_start - 1；长度 End-Start 与 SAF 的 End-Start+1 一致）。
    prom_bed <- to_bed(prom_saf)
    gb_bed   <- to_bed(gb_saf)

    write.table(prom_bed, file = file.path(argv$outdir, "promoter.bed"), sep = "\t",
                quote = FALSE, row.names = FALSE, col.names = FALSE)
    message("[gtf2saf] Promoter BED: ", nrow(prom_bed), " regions → ",
            file.path(argv$outdir, "promoter.bed"))

    write.table(gb_bed, file = file.path(argv$outdir, "genebody.bed"), sep = "\t",
                quote = FALSE, row.names = FALSE, col.names = FALSE)
    message("[gtf2saf] Genebody BED: ", nrow(gb_bed), " regions → ",
            file.path(argv$outdir, "genebody.bed"))

    # ── 分支 B：所有 transcript 的 gene body union → SAF ──
    union_saf <- build_genebody_union_saf(transcript_gr, genebody_offset)
    write.table(union_saf, file = file.path(argv$outdir, "genebody_union.saf"),
                sep = "\t", quote = FALSE, row.names = FALSE)
    message("[gtf2saf] Genebody union SAF: ", nrow(union_saf), " intervals → ",
            file.path(argv$outdir, "genebody_union.saf"))

    # ── gene 级产物（gene.bed，供 tss_meta + 信号表 intersect）──
    gene_bed <- data.frame(
        Chr = as.character(seqnames(gene_gr)), Start = start(gene_gr) - 1, End = end(gene_gr),
        GeneID = mcols(gene_gr)$gene_id, Score = ".", Strand = as.character(strand(gene_gr)),
        stringsAsFactors = FALSE
    )
    write.table(gene_bed, file = file.path(argv$outdir, "gene.bed"), sep = "\t",
                quote = FALSE, row.names = FALSE, col.names = FALSE)
    message("[gtf2saf] Gene BED: ", nrow(gene_bed), " regions → ",
            file.path(argv$outdir, "gene.bed"))

    # plus_genes.bed / minus_genes.bed 未使用，先注释掉
    # plus_idx  <- gene_bed$Strand == "+"
    # minus_idx <- gene_bed$Strand == "-"
    # write.table(gene_bed[plus_idx, ],  file = file.path(argv$outdir, "plus_genes.bed"), sep = "\t",
    #             quote = FALSE, row.names = FALSE, col.names = FALSE)
    # write.table(gene_bed[minus_idx, ], file = file.path(argv$outdir, "minus_genes.bed"), sep = "\t",
    #             quote = FALSE, row.names = FALSE, col.names = FALSE)
    # message("[gtf2saf] Strand-split BED: + ", sum(plus_idx), " genes, - ", sum(minus_idx), " genes")

    message("[gtf2saf] Done.")
}

main(argv)
