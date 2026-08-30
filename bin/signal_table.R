#!/usr/bin/env Rscript
# =============================================================================
# signal_table.R — Build the per-base Pol II active-site signal table
#
# 输入：
#   --manifest    sample \t plus_bedgraph \t minus_bedgraph（样本名来自 meta，不解析文件名）
#   --groups      YAML：group -> [samples]
#   --gene_bed    标准 BED6（Chr Start End GeneID . Strand，0-based 半开）
#   --gtf         参考 GTF（读 transcript 记录取每 gene 代表 transcript）
#   --annotation  gene 注释表（首列 gene_id，元数据列透传）
# 输出：
#   chrom start end {group}.Count {group}.Signal ... transcript_id gene_id Strand <注释列...>
#
# 约定：
#   - bedGraph 0-based 半开（genomecov -5 单碱基区间，end = start + 1），直接透传
#   - _plus = + 链基因信号（正）；_minus = - 链基因信号（负，取 abs）
#   - RPM = Count / 该组总 5' 端数 × 1e6
#   - transcript_id = 该 gene 最长 transcript（tie-breaking 与 gtf2saf.R 一致）
# =============================================================================
suppressWarnings(suppressMessages({
    library(argparser)
    library(rtracklayer)
    library(GenomicRanges)
    library(yaml)
}))

# ── Parse args ──
argv <- arg_parser("Build per-base Pol II active-site signal table")
argv <- add_argument(argv, "--manifest",   help = "sample \\t plus_bg \\t minus_bg")
argv <- add_argument(argv, "--groups",     help = "Groups YAML (group -> [samples])")
argv <- add_argument(argv, "--gene_bed",   help = "Gene BED6 (Chr Start End GeneID . Strand)")
argv <- add_argument(argv, "--gtf",        help = "Reference GTF")
argv <- add_argument(argv, "--annotation", help = "Gene annotation table (first col = gene_id)")
argv <- add_argument(argv, "--output",     help = "Output per-base signal table (tsv)")
argv <- parse_args(argv)

# ── 读 bedGraph（0-based 半开，4 列），并把 -bg 合并的区间展开为单碱基 ──
read_bedgraph <- function(f) {
    if (!file.exists(f)) stop("[signal_table] bedGraph not found: ", f)
    bg <- read.delim(f, header = FALSE, sep = "\t", stringsAsFactors = FALSE,
                     col.names = c("chrom", "start", "end", "value"))
    if (nrow(bg) == 0) return(bg)
    # genomecov -bg 会把相邻同值碱基合并成 [start,end) 区间，这里展开为单碱基（end=start+1）
    lens <- bg$end - bg$start
    data.frame(
        chrom = rep(bg$chrom, lens),
        start = rep(bg$start, lens) + sequence(lens) - 1,
        end   = rep(bg$start, lens) + sequence(lens),
        value = rep(bg$value, lens),
        stringsAsFactors = FALSE
    )
}

# ── 组内逐碱基求和（plus/minus 分别）──
aggregate_bedgraphs <- function(files) {
    all <- do.call(rbind, lapply(files, read_bedgraph))
    if (is.null(all) || nrow(all) == 0)
        return(data.frame(chrom = character(), start = integer(),
                          end = integer(), value = numeric()))
    agg <- aggregate(value ~ chrom + start + end, data = all, FUN = sum)
    agg$chrom <- as.character(agg$chrom)   # aggregate 把 chrom 转成 factor，转回字符
    agg
}

# ── 读 gene.bed（BED6）──
read_gene_bed <- function(f) {
    read.delim(f, header = FALSE, sep = "\t", stringsAsFactors = FALSE,
               col.names = c("chrom", "start", "end", "gene_id", "score", "strand"))
}

# ── 从 GTF 选每 gene 代表（最长）transcript（tie-breaking 同 gtf2saf.R）──
build_rep_transcript <- function(gtf_file) {
    gr <- rtracklayer::import(gtf_file, format = "gtf")
    gr <- gr[gr$type == "transcript"]
    gene_id <- mcols(gr)$gene_id
    tid     <- mcols(gr)$transcript_id
    strand  <- as.character(strand(gr))
    start   <- start(gr)
    end     <- end(gr)
    span    <- end - start                          # 长度量度（1-based 下真实长度 = span+1，排序等价）
    tss_key <- ifelse(strand == "+", start, -end)   # 最上游 TSS tie-break key
    if (is.null(tid) || all(is.na(tid))) {
        tid <- paste0(as.character(seqnames(gr)), ":", start, "-", end)
    } else {
        tid[is.na(tid)] <- paste0(as.character(seqnames(gr))[is.na(tid)], ":",
                                  start[is.na(tid)], "-", end[is.na(tid)])
    }
    df <- data.frame(i = seq_along(gr), gene_id = gene_id, span = span,
                     tss_key = tss_key, tid = tid, stringsAsFactors = FALSE)
    df <- df[order(df$gene_id, -df$span, df$tss_key, df$tid), ]
    keep <- !duplicated(df$gene_id)
    data.frame(gene_id = df$gene_id[keep], transcript_id = df$tid[keep],
               stringsAsFactors = FALSE)
}

# ── 单组：链向 intersect 逐碱基信号 → 每碱基 gene_id / strand / count ──
annotate_group <- function(plus_agg, minus_agg, gene_bed) {
    parts <- list()
    if (nrow(plus_agg) > 0) {
        sig <- GRanges(plus_agg$chrom, IRanges(plus_agg$start + 1, plus_agg$end))
        g   <- gene_bed[gene_bed$strand == "+", , drop = FALSE]
        ggr <- GRanges(g$chrom, IRanges(g$start + 1, g$end))
        hits <- findOverlaps(sig, ggr)
        q <- queryHits(hits); s <- subjectHits(hits)
        parts[["+"]] <- data.frame(
            chrom = plus_agg$chrom[q], start = plus_agg$start[q], end = plus_agg$end[q],
            count = plus_agg$value[q], strand = "+", gene_id = g$gene_id[s],
            stringsAsFactors = FALSE)
    }
    if (nrow(minus_agg) > 0) {
        sig <- GRanges(minus_agg$chrom, IRanges(minus_agg$start + 1, minus_agg$end))
        g   <- gene_bed[gene_bed$strand == "-", , drop = FALSE]
        ggr <- GRanges(g$chrom, IRanges(g$start + 1, g$end))
        hits <- findOverlaps(sig, ggr)
        q <- queryHits(hits); s <- subjectHits(hits)
        parts[["-"]] <- data.frame(
            chrom = minus_agg$chrom[q], start = minus_agg$start[q], end = minus_agg$end[q],
            count = abs(minus_agg$value[q]), strand = "-", gene_id = g$gene_id[s],
            stringsAsFactors = FALSE)
    }
    do.call(rbind, parts)
}

# ── Main ──
main <- function(argv) {
    manifest <- read.delim(argv$manifest, header = FALSE, sep = "\t",
                           stringsAsFactors = FALSE, col.names = c("sample", "plus", "minus"))
    groups   <- yaml::read_yaml(argv$groups)
    gene_bed <- read_gene_bed(argv$gene_bed)
    rep_tr   <- build_rep_transcript(argv$gtf)

    # ── 每组：逐碱基聚合 + 链向注释 ──
    group_tables <- list()
    group_totals <- c()
    for (grp in names(groups)) {
        samples <- as.character(unlist(groups[[grp]]))
        m <- manifest[manifest$sample %in% samples, , drop = FALSE]
        if (nrow(m) == 0)
            stop("[signal_table] group '", grp, "' has no samples in manifest")

        plus_agg  <- aggregate_bedgraphs(m$plus)
        minus_agg <- aggregate_bedgraphs(m$minus)

        total <- sum(plus_agg$value) + sum(abs(minus_agg$value))   # 组总 5' 端数
        group_totals[grp] <- total

        gt <- annotate_group(plus_agg, minus_agg, gene_bed)
        if (is.null(gt) || nrow(gt) == 0)
            gt <- data.frame(chrom = character(), start = integer(), end = integer(),
                             count = numeric(), strand = character(), gene_id = character())
        gt <- gt[!duplicated(gt), ]
        names(gt)[names(gt) == "count"] <- paste0(grp, ".Count")
        group_tables[[grp]] <- gt
    }

    # ── 跨组合并（full outer join on 碱基+基因），缺失组填 0 ──
    key_cols <- c("chrom", "start", "end", "strand", "gene_id")
    all_keys <- unique(do.call(rbind, lapply(group_tables, function(d) d[, key_cols, drop = FALSE])))
    if (is.null(all_keys) || nrow(all_keys) == 0) {
        write.table(data.frame(), file = argv$output, sep = "\t", quote = FALSE,
                    row.names = FALSE, col.names = FALSE)
        message("[signal_table] no signal -> empty output")
        return(invisible(NULL))
    }
    merged <- all_keys
    for (grp in names(group_tables)) {
        gt <- group_tables[[grp]]
        cnt_col <- paste0(grp, ".Count")
        merged <- merge(merged, gt[, c(key_cols, cnt_col), drop = FALSE],
                        by = key_cols, all.x = TRUE, sort = FALSE)
        merged[[cnt_col]][is.na(merged[[cnt_col]])] <- 0
        sig_col <- paste0(grp, ".Signal")
        if (group_totals[grp] > 0) {
            merged[[sig_col]] <- merged[[cnt_col]] / group_totals[grp] * 1e6
        } else {
            merged[[sig_col]] <- 0   # 组总 5' 端数为 0（无信号）时 Signal 记 0，避免 0/0
        }
    }

    # ── 拼 transcript_id + gene 注释 ──
    merged <- merge(merged, rep_tr, by = "gene_id", all.x = TRUE, sort = FALSE)
    merged$transcript_id[is.na(merged$transcript_id)] <- ""

    ann <- read.delim(argv$annotation, header = TRUE, sep = "\t",
                      stringsAsFactors = FALSE, check.names = FALSE,
                      quote = "", comment.char = "")
    ann[is.na(ann)] <- ""
    names(ann)[1] <- "gene_id"
    merged <- merge(merged, ann, by = "gene_id", all.x = TRUE, sort = FALSE)

    # ── 列序：chrom start end {group}.Count {group}.Signal ... transcript_id gene_id Strand <注释列> ──
    names(merged)[names(merged) == "strand"] <- "Strand"
    group_cols <- unlist(lapply(names(group_tables),
                                function(g) c(paste0(g, ".Count"), paste0(g, ".Signal"))))
    ann_cols <- setdiff(names(ann), "gene_id")
    out <- merged[c("chrom", "start", "end", group_cols,
                    "transcript_id", "gene_id", "Strand", ann_cols)]

    write.table(out, file = argv$output, sep = "\t", quote = FALSE, row.names = FALSE)
    message("[signal_table] ", nrow(out), " rows x ", ncol(out), " cols -> ", argv$output)
}

main(argv)
