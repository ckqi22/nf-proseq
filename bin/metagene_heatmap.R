#!/usr/bin/env Rscript
# =============================================================================
# metagene_heatmap.R — Per-gene TSS-centered heatmap from strand-specific PRO-seq
# bigWigs. Each row is a gene, sorted by total window signal (strongest at the
# bottom); columns are bins from --upstream to --downstream relative to the TSS.
#
# 以每个基因 TSS 为锚点提取 [--upstream, --downstream] 窗口的 sense 链信号，
# 每行一个基因，按窗口内总信号升序排序（信号最强者排在底部，comet 样式）。
#
# 输入约定（与 metagene_profile.R 一致）：
#   - _plus  = + 链基因信号（正值）；_minus = - 链基因信号（负值，读取后 abs()）。
#   - bigWig 已 CPM 归一化；本脚本按 RPKM(= CPM × 1000) 展示，保留原 dheatmap 色阶。
#
# Input:
#   --plus_bw             + strand bigWig (CPM). Required.
#   --minus_bw            - strand bigWig (CPM, negative values). Required.
#   --tss_bed             TSS BED6（每行一个待绘图 gene/feature；end=TSS 1-based）. Required.
#   --outdir              输出目录，默认 ./
#   --prefix              输出文件名前缀（区分样本），默认无
#   --upstream            TSS 上游 bp，默认 50
#   --downstream          TSS 下游 bp，默认 150
#   --bin                 分箱大小 bp，默认 1（单碱基）
#   --max_cut             色阶上限 = 数据最大值 × max_cut，默认 1.0
#   --smooth              滑动平均平滑窗口 bp（奇数；0 = 不平滑），默认 0
#   --min_signal          仅保留窗口内总 RPKM 不低于该值的基因（去掉无信号行），默认 0
#
# Output:
#   <prefix>_metagene_heatmap_matrix.tsv   gene_id / <每 bin 信号(RPKM)>
#                                          （平滑后、已过滤/排序，行序与热图一致）
#   <prefix>_metagene_heatmap.png
#
# Usage:
#   Rscript metagene_heatmap.R \
#       --plus_bw sample_plus.bigWig --minus_bw sample_minus.bigWig \
#       --tss_bed tss.bed --prefix S1
# =============================================================================
suppressWarnings(suppressMessages({
    library(argparser)
    library(GenomicRanges)
    library(IRanges)
    library(GenomeInfoDb)
    library(rtracklayer)
    library(ggplot2)
}))

# ── Parse args ──
argv <- arg_parser("Per-gene TSS heatmap from strand-specific PRO-seq bigWigs")
argv <- add_argument(argv, "--plus_bw", help = "+ strand bigWig (CPM). Required.")
argv <- add_argument(argv, "--minus_bw", help = "- strand bigWig (CPM, negative). Required.")
argv <- add_argument(argv, "--tss_bed", help = "TSS BED6 (one line per feature to plot). Required.")
argv <- add_argument(argv, "--outdir", help = "Output directory", default = "./")
argv <- add_argument(argv, "--prefix", help = "Output filename prefix to distinguish samples", default = "")
argv <- add_argument(argv, "--upstream", help = "Bases upstream of TSS", default = 50, type = "integer")
argv <- add_argument(argv, "--downstream", help = "Bases downstream of TSS", default = 150, type = "integer")
argv <- add_argument(argv, "--bin", help = "Bin size in bp", default = 1, type = "integer")
argv <- add_argument(argv, "--max_cut", help = "Color upper limit = data max x max_cut", default = 1.0, type = "double")
argv <- add_argument(argv, "--smooth", help = "Rolling-mean smoothing window in bp (odd; 0 = off)", default = 101, type = "integer")
argv <- add_argument(argv, "--min_signal", help = "Drop genes with total window RPKM below this", default = 0, type = "double")
argv <- parse_args(argv)

# —— test args ——
# argv$plus_bw <- "/workplace/chenkai/proseq_test/GSE181161/results/08.Pol2_coverage/HS0_rep1_plus_cpm.bigWig"
# argv$minus_bw <- "/workplace/chenkai/proseq_test/GSE181161/results/08.Pol2_coverage/HS0_rep1_minus_cpm.bigWig"
# argv$tss_bed <- "path/to/tss.bed"   # 由 gtf2bed.R 生成

# ── Read TSS from a BED6 (one line per feature; end = 1-based TSS) ──
read_tss_bed <- function(f) {
    message("[metagene_heatmap] Reading TSS BED: ", f)
    bed <- read.delim(f, header = FALSE, sep = "\t", comment.char = "#",
                      stringsAsFactors = FALSE,
                      col.names = c("chrom", "start", "end", "name", "score", "strand"))
    if (nrow(bed) == 0)
        stop("[metagene_heatmap] TSS BED is empty: ", f)
    if (!all(bed$strand %in% c("+", "-")))
        stop("[metagene_heatmap] TSS BED strand column must be '+' or '-'")
    tss <- bed$end                          # 1-based TSS（BED 半开 [start,end)，end=TSS）
    GRanges(seqnames = bed$chrom,
            ranges   = IRanges(start = tss, end = tss),
            strand   = bed$strand,
            gene_id  = bed$name)
}

# ── Bin a gene × position matrix (columns) into bins of width `bin` (mean per bin) ──
bin_matrix <- function(mat, bin) {
    if (bin <= 1) return(mat)
    n    <- ncol(mat)
    nbin <- n / bin                       # (--upstream + --downstream) %% --bin == 0 已在 main 校验
    out  <- matrix(NA_real_, nrow = nrow(mat), ncol = nbin)
    for (j in seq_len(nbin)) {
        cc <- ((j - 1) * bin + 1):(j * bin)
        out[, j] <- rowMeans(mat[, cc, drop = FALSE])
    }
    out
}

# ── 滑动平均平滑（按行，边缘截断）；w 为 bin 数，奇数居中对称。w<=1 不平滑 ──
smooth_matrix <- function(mat, w) {
    if (is.null(w) || w <= 1) return(mat)
    if (w %% 2 == 0) w <- w + 1                 # 保证奇数，窗口以当前 bin 为中心
    half <- w %/% 2
    n    <- ncol(mat)
    if (w >= n) return(matrix(rowMeans(mat, na.rm = TRUE),
                              nrow = nrow(mat), ncol = n))
    out <- mat
    for (i in seq_len(n)) {
        lo <- max(1L, i - half)
        hi <- min(n, i + half)
        out[, i] <- rowMeans(mat[, lo:hi, drop = FALSE], na.rm = TRUE)
    }
    out
}

# white rainbow 色带（与原 dheatmap 的 heatmap.h mode 9 的 9 个节点色一致）
heatmap_colors <- c("#FFFFFF", "#EE82EE", "#4B0082", "#0000FF", "#005F80",
                    "#008000", "#FFFF00", "#FFA500", "#FF0000")

# ── Main ──
main <- function(argv) {
    # ── validate inputs ──
    if (is.null(argv$plus_bw) || is.null(argv$minus_bw))
        stop("[metagene_heatmap] both --plus_bw and --minus_bw are required")
    if (is.null(argv$tss_bed))
        stop("[metagene_heatmap] --tss_bed is required")
    if (!file.exists(argv$plus_bw))
        stop("[metagene_heatmap] --plus_bw not found: ", argv$plus_bw)
    if (!file.exists(argv$minus_bw))
        stop("[metagene_heatmap] --minus_bw not found: ", argv$minus_bw)
    if (!file.exists(argv$tss_bed))
        stop("[metagene_heatmap] --tss_bed not found: ", argv$tss_bed)

    U <- as.integer(argv$upstream)
    D <- as.integer(argv$downstream)
    B <- as.integer(argv$bin)
    if (U <= 0) stop("[metagene_heatmap] --upstream must be > 0")
    if (D <= 0) stop("[metagene_heatmap] --downstream must be > 0")
    if (B <= 0) stop("[metagene_heatmap] --bin must be > 0")
    if (!is.finite(argv$max_cut) || argv$max_cut <= 0)
        stop("[metagene_heatmap] --max_cut must be > 0")
    if (!is.finite(argv$smooth) || argv$smooth < 0)
        stop("[metagene_heatmap] --smooth must be >= 0")
    if (!is.finite(argv$min_signal) || argv$min_signal < 0)
        stop("[metagene_heatmap] --min_signal must be >= 0")
    if ((U + D) %% B != 0)
        stop("[metagene_heatmap] (--upstream + --downstream) must be divisible by --bin")

    W  <- U + D                          # 窗口总宽（bp）
    NB <- W / B                          # 分箱数

    # ── build TSS GRanges ──
    tss_gr <- read_tss_bed(argv$tss_bed)

    # seqlengths 取自 bigWig（BED 可能缺失 / 含 scaffold）
    bw_si  <- seqinfo(rtracklayer::BigWigFile(argv$plus_bw))
    common <- intersect(seqlevels(tss_gr), seqlevels(bw_si))
    if (length(common) == 0)
        stop("[metagene_heatmap] no shared seqlevels between BED and bigWig; check chromosome naming")
    tss_gr <- keepSeqlevels(tss_gr, common, pruning.mode = "coarse")
    seqlengths(tss_gr) <- seqlengths(bw_si)[common]
    message("[metagene_heatmap] Using ", length(tss_gr), " representative TSS")

    plus_bw  <- rtracklayer::BigWigFile(argv$plus_bw)
    minus_bw <- rtracklayer::BigWigFile(argv$minus_bw)

    # ── per-gene window extraction（5'->3' 定向）──
    n_genes <- length(tss_gr)
    mat     <- matrix(NA_real_, nrow = n_genes, ncol = W)
    used    <- logical(n_genes)

    chroms <- as.character(unique(seqnames(tss_gr)))
    for (ch in chroms) {
        idx <- which(as.character(seqnames(tss_gr)) == ch)
        if (length(idx) == 0) next
        sl <- seqlengths(tss_gr)[ch]
        if (is.na(sl)) sl <- seqlengths(seqinfo(plus_bw))[ch]   # 兜底：bigWig 染色体长度
        which_gr <- GRanges(ch, IRanges(1, sl))
        plus_rle  <- rtracklayer::import(plus_bw,  which = which_gr, as = "RleList")[[ch]]
        minus_rle <- abs(rtracklayer::import(minus_bw, which = which_gr, as = "RleList")[[ch]])
        for (i in idx) {
            gstr <- as.character(strand(tss_gr)[i])
            pos  <- start(tss_gr)[i]
            if (gstr == "+") {
                s <- pos - U; e <- pos + D - 1L          # 5'->3' 与基因组方向一致
                if (s < 1 || e > sl) next
                v <- as.numeric(plus_rle[s:e])
            } else {
                s <- pos - D + 1L; e <- pos + U          # 反向取 minus 链，再 rev 回 5'->3'
                if (s < 1 || e > sl) next
                v <- rev(as.numeric(minus_rle[s:e]))
            }
            mat[i, ] <- v
            used[i]  <- TRUE
        }
    }

    n_used <- sum(used)
    message("[metagene_heatmap] Extracted ", n_used, " / ", n_genes,
            " TSS windows (others fell partly outside chromosome bounds)")
    if (n_used == 0)
        stop("[metagene_heatmap] no TSS windows could be extracted; check chromosome naming / seqlengths")

    mat       <- mat[used, , drop = FALSE]
    genes_out <- tss_gr[used]

    # ── 分箱（bin > 1 时取 bin 内 mean）──
    mat <- bin_matrix(mat, B)

    # ── CPM -> RPKM（×1000，与原 dheatmap.cpp 一致）──
    rpkm <- mat * 1000

    # ── 滑动平均平滑（缓解单碱基稀疏；0 = 不平滑）──
    if (argv$smooth > 0) {
        w_smooth <- as.integer(round(argv$smooth / B))   # bp -> bin 数
        rpkm <- smooth_matrix(rpkm, w_smooth)
    }

    # ── 过滤：去掉窗口内总信号过低的基因（默认 0 = 全部保留）──
    tot  <- rowSums(rpkm, na.rm = TRUE)
    keep <- tot >= argv$min_signal
    if (!any(keep))
        stop("[metagene_heatmap] no genes remain after --min_signal filter; lower --min_signal")
    if (argv$min_signal > 0)
        message("[metagene_heatmap] Kept ", sum(keep), " / ", length(keep),
                " genes after --min_signal filter")
    rpkm      <- rpkm[keep, , drop = FALSE]
    genes_out <- genes_out[keep]
    tot       <- tot[keep]

    # ── 排序：按窗口总信号升序（最弱在前，绘图时 y 反向使最强在底部）──
    ord       <- order(tot)
    rpkm      <- rpkm[ord, , drop = FALSE]
    genes_out <- genes_out[ord]

    # ── log10(RPKM+1)，再按 max_cut 封顶 ──
    mat_log <- log10(rpkm + 1)
    cap     <- max(mat_log, na.rm = TRUE) * argv$max_cut
    if (!is.finite(cap) || cap <= 0) cap <- 1             # 全零信号时避免退化色阶
    mat_log[mat_log > cap] <- cap

    # 色阶刻度标回 RPKM（0, 1, 10, 100, 1000）
    rk   <- c(0, 1, 10, 100, 1000)
    rk   <- rk[log10(rk + 1) <= cap]
    brks <- log10(rk + 1)
    labs <- parse(text = ifelse(rk > 10, paste0("10^", log10(rk)), as.character(rk)))

    # 位置（bin 左边界，相对 TSS）
    pos <- seq(-U, D - B, by = B)

    prefix <- if (is.null(argv$prefix) || argv$prefix == "") "" else paste0(argv$prefix, "_")

    # ── 写矩阵（RPKM，行序与热图一致）──
    mat_df <- data.frame(gene_id = genes_out$gene_id, stringsAsFactors = FALSE)
    mat_df[as.character(pos)] <- rpkm
    mat_out <- file.path(argv$outdir, paste0(prefix, "metagene_heatmap_matrix.tsv"))
    write.table(mat_df, file = mat_out, sep = "\t", row.names = FALSE, quote = FALSE)
    message("[metagene_heatmap] Wrote matrix: ", mat_out)

    # ── 绘图 ──
    message("[metagene_heatmap] Plotting ...")
    plot_df <- data.frame(
        gene     = rep(seq_len(nrow(mat_log)), ncol(mat_log)),
        position = rep(pos + B/2, each = nrow(mat_log)),   # bin 中心（geom_raster 单元格中心）
        signal   = as.vector(mat_log)
    )

    p <- ggplot(plot_df, aes(x = position, y = gene, fill = signal)) +
        geom_raster(interpolate = FALSE) +
        scale_fill_gradientn(colors = heatmap_colors, limits = c(0, cap),
                             breaks = brks, labels = labs,
                             name = "PRO-seq (RPKM)",
                             guide = guide_colorbar(barheight = unit(40, "mm"),
                                                    barwidth  = unit(3, "mm"),
                                                    ticks     = FALSE)) +
        scale_x_continuous(expand = c(0, 0)) +
        scale_y_reverse(expand = c(0, 0)) +
        geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
        labs(x = "Position from TSS (bp)", y = NULL) +
        theme_classic(base_family = "Liberation Sans", base_size = 12) +
        theme(axis.line.y = element_blank(),
              axis.text.y  = element_blank(),
              axis.ticks.y = element_blank(),
              legend.position = "right")

    png_out <- file.path(argv$outdir, paste0(prefix, "metagene_heatmap.png"))
    ggsave(png_out, p, width = 5.2, height = 4.5, dpi = 300, units = "in")
    message("[metagene_heatmap] Wrote figure: ", png_out)
    message("[metagene_heatmap] Done.")
}

main(argv)
