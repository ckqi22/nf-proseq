#!/usr/bin/env Rscript
# =============================================================================
# spikein_scale.R — 由 spike-in read 数计算每样本缩放因子
#
# 公式：
#   factor_i      = 1e6 / spike_count_i            # RPM-spike：乘到 bigWig / 矩阵
#   size_factor_i = geomean(spike_count) / spike_count_i   # 相对因子（geomean=1，仅报告）
#
# 输入：
#   --counts      两列（无表头）：sample \t spike_count
#
# 输出：
#   spikein_scale_factors.tsv   sample spike_count factor size_factor
#
# Usage:
#   Rscript spikein_scale.R --counts spike_counts.tsv --output_dir ./
# =============================================================================
suppressWarnings(suppressMessages({
    library(argparser)
}))

argv <- arg_parser("Spike-in scale factors + matrix scaling")
argv <- add_argument(argv, "--counts",     help = "Two columns (no header): sample \\t spike_count")
argv <- add_argument(argv, "--output_dir", help = "Output directory", default = "./")
argv <- parse_args(argv)

main <- function(argv) {
    if (is.null(argv$counts) || !file.exists(argv$counts))
        stop("[spikein_scale] --counts is required")

    ct <- read.delim(argv$counts, header = FALSE, sep = "\t", stringsAsFactors = FALSE,
                     col.names = c("sample", "spike_count"))
    dir.create(argv$output_dir, showWarnings = FALSE, recursive = TRUE)

    # spike_count ≤ 0 或 NA：直接报错（不写 Inf/NA 继续；上游 spikein.nf 已按 fraction 拦截，此处兜底）
    bad <- is.na(ct$spike_count) | ct$spike_count <= 0
    if (any(bad))
        stop("[spikein_scale] spike_count <= 0 or NA for sample(s): ",
             paste(ct$sample[bad], collapse = ", "))

    geomean <- exp(mean(log(ct$spike_count)))

    ct$factor      <- 1e6 / ct$spike_count
    ct$size_factor <- geomean / ct$spike_count

    out_factors <- file.path(argv$output_dir, "spikein_scale_factors.tsv")
    write.table(ct, out_factors, sep = "\t", row.names = FALSE, quote = FALSE)
    message("[spikein_scale] wrote ", out_factors)
}

main(argv)
