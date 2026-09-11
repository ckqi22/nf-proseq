#!/usr/bin/env Rscript
# ============================================================
# overall_statistics.R — 合并 read_statistics + bowtie2 比对率
#
# 参照 TT-seq 的 overall_statistics.xlsx（read stats + mapping stats 合并表）。
# nf-proseq 没有单独的 mapping_statistics.txt，只有 bowtie2 的 stderr 摘要
# （ALIGN 里 `2> >(tee ${sample}_summary_bowtie2.txt)`），这里直接从中取
# 「overall alignment rate」合并进 overall_statistics。
#
# 用法:
#   Rscript overall_statistics.R \
#       --read_statistics     results/03.Data_QC/read_statistics.txt \
#       --bowtie2_summary_dir results/04.Alignment \
#       --output              results/03.Data_QC/overall_statistics.txt
#
# 输出列 = Sample + read_statistics 全部列 + Overall Mapped Rate（0~1 小数，
#          与 Clean Ratio 一致，后续 txt2xlsx 会按 Rate 列渲染成百分比）。
# ============================================================

# ---- 极简参数解析（无第三方依赖）----
argv <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default = NULL) {
  idx <- which(argv == name)
  if (length(idx) == 0) return(default)
  if (idx + 1 > length(argv)) stop("缺少参数值: ", name)
  argv[idx + 1]
}
read_statistics <- get_arg("--read_statistics")
summary_dir     <- get_arg("--bowtie2_summary_dir")
output          <- get_arg("--output", "overall_statistics.txt")

if (is.null(read_statistics) || is.null(summary_dir))
  stop("需要 --read_statistics 与 --bowtie2_summary_dir")

if (!file.exists(read_statistics)) stop("read_statistics.txt 不存在: ", read_statistics)

# ---- 读取 read_statistics ----
stat <- read.table(read_statistics, header = TRUE, sep = "\t", quote = "",
                   check.names = FALSE, stringsAsFactors = FALSE, fill = TRUE)
if (!"Sample" %in% colnames(stat)) stop("read_statistics 缺 Sample 列")

# ---- 从 bowtie2 摘要取 overall alignment rate ----
summary_files <- list.files(summary_dir, pattern = "_summary_bowtie2\\.txt$",
                            full.names = TRUE)
if (length(summary_files) == 0)
  stop("bowtie2 摘要目录下无 *_summary_bowtie2.txt: ", summary_dir)

extract_rate <- function(f) {
  lines <- readLines(f, warn = FALSE)
  hits  <- grep("overall alignment rate", lines, value = TRUE)
  if (length(hits) == 0) return(NA_real_)
  pct <- sub("^[^0-9.]*([0-9.]+)%.*$", "\\1", hits[1])
  as.numeric(pct) / 100
}

rate_map <- setNames(rep(NA_real_, length(summary_files)),
                     sub("_summary_bowtie2\\.txt$", "", basename(summary_files)))
for (f in summary_files) {
  smp <- sub("_summary_bowtie2\\.txt$", "", basename(f))
  rate_map[[smp]] <- extract_rate(f)
}

# 匹配 read_statistics 的样本顺序；缺失样本填 NA 并告警
mapped <- rate_map[stat$Sample]
if (anyNA(mapped)) {
  miss <- stat$Sample[is.na(mapped)]
  message("[overall_statistics] 以下样本缺 bowtie2 摘要，Overall Mapped Rate 置 NA: ",
          paste(miss, collapse = ", "))
}

out <- cbind(stat, `Overall Mapped Rate` = unname(mapped))
write.table(out, file = output, sep = "\t", row.names = FALSE, quote = FALSE)
message("[overall_statistics] 写出: ", output, " (", nrow(out), " 样本 x ",
        ncol(out), " 列)")
