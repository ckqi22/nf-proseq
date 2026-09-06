#!/usr/bin/env Rscript
# =============================================================================
# metagene_profile.R — TSS-centered metagene profile (±WINDOW bp) from
# strand-specific PRO-seq bigWigs, showing sense / antisense signal on one plot.
#
# 将每个基因按自身链取向到 5'->3'，绘制 TSS 附近的 sense（正义）/ antisense（反义）信号。
# 负链存负值，读取后 abs()。
#
# Input:
#   --plus_bw    + strand bigWig. Required.
#   --minus_bw   - strand bigWig. Required.
#   --tss_bed    TSS BED6（每行一个待绘图 gene/feature；end=TSS 1-based）. Required.
#   --prefix     输出文件名前缀（区分样本），默认无前缀
#
# Output:
#   <prefix>_metagene_profile_matrix.tsv  position / sense / antisense per bin
#   <prefix>_metagene_profile.png
#
# Usage:
#   Rscript metagene_profile.R \
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
argv <- arg_parser("TSS-centered metagene profile from strand-specific PRO-seq bigWigs")
argv <- add_argument(argv, "--plus_bw", help = "+ strand bigWig. Required.")
argv <- add_argument(argv, "--minus_bw", help = "- strand bigWig. Required.")
argv <- add_argument(argv, "--tss_bed", help = "TSS BED6 (one line per feature to plot). Required.")
argv <- add_argument(argv, "--outdir", help = "Output directory", default = "./")
argv <- add_argument(argv, "--prefix", help = "Output filename prefix to distinguish samples", default = "")
argv <- add_argument(argv, "--window", help = "Half-window around TSS in bp", default = 1000, type = "integer")
argv <- add_argument(argv, "--bin", help = "Bin size in bp", default = 10, type = "integer")
argv <- add_argument(argv, "--mirror", flag = TRUE, help = "Mirror style: antisense plotted as negative values below zero")
argv <- parse_args(argv)

# ── Read TSS from a BED6 (one line per feature; end = 1-based TSS) ──
read_tss_bed <- function(f) {
    message("[metagene_profile] Reading TSS BED: ", f)
    bed <- read.delim(f, header = FALSE, sep = "\t", comment.char = "#",
                      stringsAsFactors = FALSE,
                      col.names = c("chrom", "start", "end", "name", "score", "strand"))
    if (nrow(bed) == 0)
        stop("[metagene_profile] TSS BED is empty: ", f)
    if (!all(bed$strand %in% c("+", "-")))
        stop("[metagene_profile] TSS BED strand column must be '+' or '-'")
    tss <- bed$end                          # 1-based TSS（BED 半开 [start,end)，end=TSS）
    GRanges(seqnames = bed$chrom,
            ranges   = IRanges(start = tss, end = tss),
            strand   = bed$strand)
}

# ── Bin a length-2W numeric vector into NB bins of width B (mean per bin) ──
binvec <- function(v, B) colMeans(matrix(v, nrow = B))

# ── Main ──
main <- function(argv) {
    # ── validate inputs ──
    if (is.null(argv$plus_bw) || is.null(argv$minus_bw))
        stop("[metagene_profile] both --plus_bw and --minus_bw are required")
    if (is.null(argv$tss_bed))
        stop("[metagene_profile] --tss_bed is required")

    W   <- as.integer(argv$window)   # half window
    B   <- as.integer(argv$bin)      # bin size
    NB  <- as.integer(2 * W / B)     # number of bins (2W must be divisible by B)
    if ((2 * W) %% B != 0) stop("[metagene_profile] 2*window must be divisible by --bin")
    stopifnot(NB > 0)

    # ── build TSS GRanges ──
    if (!file.exists(argv$tss_bed))
        stop("[metagene_profile] --tss_bed not found: ", argv$tss_bed)
    tss_gr <- read_tss_bed(argv$tss_bed)

    # set seqlengths from the bigWig (BED may lack them / have scaffolds)
    bw_si  <- seqinfo(rtracklayer::BigWigFile(argv$plus_bw))
    common <- intersect(seqlevels(tss_gr), seqlevels(bw_si))
    tss_gr <- keepSeqlevels(tss_gr, common, pruning.mode = "coarse")
    seqlengths(tss_gr) <- seqlengths(bw_si)[common]
    message("[metagene_profile] Using ", length(tss_gr), " representative TSS")

    # ── per-bin sense / antisense profiles ──
    # deepTools bamCoverage stores reverse strand as NEGATIVE values -> use abs()
    plus_bw  <- rtracklayer::BigWigFile(argv$plus_bw)
    minus_bw <- rtracklayer::BigWigFile(argv$minus_bw)

    sense_sum     <- numeric(NB)
    antisense_sum <- numeric(NB)
    n_used <- 0L

    chroms <- as.character(unique(seqnames(tss_gr)))
    for (ch in chroms) {
        idx <- which(as.character(seqnames(tss_gr)) == ch)
        if (length(idx) == 0) next
        # import full-chromosome Rle for this chrom (PRO-seq is sparse -> cheap)
        sl <- seqlengths(tss_gr)[ch]
        if (is.na(sl)) {
            # fall back to bigWig chrom length
            sl <- seqlengths(seqinfo(rtracklayer::BigWigFile(argv$plus_bw)))[ch]
        }
        which_gr <- GRanges(ch, IRanges(1, sl))
        plus_rle  <- rtracklayer::import(plus_bw,  which = which_gr, as = "RleList")[[ch]]
        minus_rle <- rtracklayer::import(minus_bw, which = which_gr, as = "RleList")[[ch]]
        minus_rle <- abs(minus_rle)   # reverse strand stored negative
        for (i in idx) {
            gstr <- as.character(strand(tss_gr)[i])
            pos  <- start(tss_gr)[i]
            if (gstr == "+") {
                s <- pos - W; e <- pos + W - 1L          # width 2W, TSS at index W+1
                if (s < 1 || e > sl) next
                vp <- as.numeric(plus_rle[s:e])
                vm <- as.numeric(minus_rle[s:e])
                sense <- binvec(vp, B); antisense <- binvec(vm, B)
            } else {
                s <- pos - W + 1L; e <- pos + W          # genomic; reverse to orient 5'->3'
                if (s < 1 || e > sl) next
                vp <- as.numeric(plus_rle[s:e])
                vm <- as.numeric(minus_rle[s:e])
                sense <- binvec(rev(vm), B)              # gene strand = minus
                antisense <- binvec(rev(vp), B)          # opposite strand = plus
            }
            sense_sum <- sense_sum + sense
            antisense_sum <- antisense_sum + antisense
            n_used <- n_used + 1L
        }
    }
    message("[metagene_profile] Profiled ", n_used, " / ", length(tss_gr),
            " TSS (others fell partly outside chromosome bounds)")
    if (n_used == 0) stop("[metagene_profile] no TSS windows could be extracted; check chromosome naming / seqlengths")

    sense_mean     <- sense_sum / n_used
    antisense_mean <- antisense_sum / n_used
    if (isTRUE(argv$mirror)) antisense_mean <- -antisense_mean   # mirror: antisense below zero

    # bin midpoints relative to TSS (bp)
    mids <- seq(-W + B/2, W - B/2, by = B)
    prof <- data.frame(position = mids, sense = sense_mean, antisense = antisense_mean)

    # output filename prefix (distinguish samples); empty -> no prefix
    prefix <- if (is.null(argv$prefix) || argv$prefix == "") "" else paste0(argv$prefix, "_")

    # save matrix
    mat_out <- file.path(argv$outdir, paste0(prefix, "metagene_profile_matrix.tsv"))
    write.table(prof, file = mat_out, sep = "\t", row.names = FALSE, quote = FALSE)
    message("[metagene_profile] Wrote matrix: ", mat_out)

    # ── plot ──
    df <- data.frame(
        position = rep(prof$position, 2),
        signal   = c(prof$sense, prof$antisense),
        strand   = factor(rep(c("sense", "antisense"), each = nrow(prof)),
                          levels = c("sense", "antisense"))
    )

    p <- ggplot(df, aes(position, signal, colour = strand)) +
        geom_line(linewidth = 0.8) +
        geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40")
    if (isTRUE(argv$mirror)) {
        ylim <- max(abs(df$signal), na.rm = TRUE)
        p <- p + geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.4) +
            scale_y_continuous(limits = c(-ylim, ylim))
    }
    p <- p +
        scale_x_continuous(breaks = seq(-W, W, by = W/2),
                           labels = function(x) sprintf("%+.0f", x)) +
        labs(x = "Distance from TSS (bp)", y = "Mean signal per bin", colour = NULL) +
        theme_classic(base_family = "Liberation Sans", base_size = 12) +
        theme(panel.grid.major.y = element_line(colour = "grey92"))

    png_out <- file.path(argv$outdir, paste0(prefix, "metagene_profile.png"))
    ggsave(png_out, p, width = 6.5, height = 4.2, dpi = 300, units = "in")
    message("[metagene_profile] Wrote figure: ", png_out)
    message("[metagene_profile] Done.")
}

main(argv)
