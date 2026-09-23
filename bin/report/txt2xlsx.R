#!/usr/bin/env Rscript
# ============================================================
# TT-Seq 通用 txt -> xlsx 转换 (openxlsx) —— 唯一 xlsx 风格源
# 交付格式: 标题 + 注释 + 彩色表头（仿公司样式）+ 数字格式
#   - 简单表 (read_statistics 等): 表头统一浅蓝 #ADD8E6
#   - 表达谱 (.Count/.Corrected/.Fpkm/.Cpm/.Signal 分组): 多色表头
#   - --group_by <列>: 按该列分 sheet（每个取值一个 sheet，该列不进 sheet）
#   - --sort_by <列>: 按该列排序（如 Sample）
#   - 自动数字格式：列名匹配 Ratio|Rate|Q30|percentage → 百分比；Reads|Bases → 千分位
#   - --diff_dir + --diff_suffix: 差异表多文件合并模式
#       （glob 目录下 *{diff_suffix} 文件，每个按 --group_by(Regulation) 拆 up/down sheet，
#         sheet 名 = "{up|down}.{文件名去后缀}"，再用 --sheet_strip 去掉尾部 _ 段）
# ============================================================
suppressPackageStartupMessages({
  library(optparse)
  library(openxlsx)
})

option_list <- list(
  make_option(c("--input"),  type = "character", help = "输入 txt (tab 分隔, 带表头)"),
  make_option(c("--output"), type = "character", help = "输出 xlsx 路径"),
  make_option(c("--title"),  type = "character", default = "", help = "标题（置顶合并加粗）"),
  make_option(c("--note"),   type = "character", default = "", help = "注释（标题下方合并自动换行）"),
  make_option(c("--note_file"), type = "character", default = "", help = "从文件读取注释（--note 为空时生效）"),
  make_option(c("--title_align"), type = "character", default = "left", help = "标题对齐方式 left/center/right（默认 left）"),
  make_option(c("--group_by"), type = "character", default = "", help = "按此列分 sheet（每个取值一个 sheet）"),
  make_option(c("--sort_by"),  type = "character", default = "", help = "按此列排序（如 Sample）"),
  make_option(c("--diff_dir"),    type = "character", default = "", help = "差异表多文件合并模式：输入目录"),
  make_option(c("--diff_suffix"), type = "character", default = "", help = "差异表文件后缀（如 _all.txt / _de.txt）"),
  make_option(c("--sheet_strip"), type = "integer", default = 0, help = "sheet 名去掉文件名尾部 N 个 _ 段（差异表用）")
)
opt <- parse_args(OptionParser(option_list = option_list))

if (!nzchar(opt$note) && nzchar(opt$note_file) && file.exists(opt$note_file)) {
  opt$note <- paste(readLines(opt$note_file, warn = FALSE), collapse = "\n")
}

read_df <- function(path) {
  read.table(path, header = TRUE, sep = "\t", quote = "", comment.char = "",
             check.names = FALSE, stringsAsFactors = FALSE)
}

# 表头颜色：按列名后缀分组，无后缀的落 default（浅蓝）
group_of <- function(name) {
  if (grepl("\\.Corrected$", name)) "corrected"
  else if (grepl("\\.Count$",  name)) "count"
  else if (grepl("\\.Fpkm$",    name)) "fpkm"
  else if (grepl("\\.Signal$",  name)) "signal"
  else if (grepl("\\.Cpm$",     name)) "cpm"
  else "default"
}
palette <- c(default  = "#ADD8E6",  # 浅蓝（Gene_id / 注释列 / 简单表）
            count     = "#FFFF00",  # 黄
            corrected = "#F4B183",  # 橙
            fpkm      = "#A9D08E",  # 绿
            signal    = "#B39DDB",  # 紫（.Signal 独立新色）
            cpm       = "#D9D9D9")  # 灰

wb <- createWorkbook()
modifyBaseFont(wb, fontSize = 11, fontName = "Times New Roman")

# ---------- 写单个 sheet ----------
write_sheet <- function(wb, sheet_name, sub, title, note) {
  nc <- ncol(sub)
  addWorksheet(wb, sheet_name)
  cur_row <- 1

  # --- 标题 ---
  if (nzchar(title)) {
    writeData(wb, sheet_name, data.frame(x = title),
              startCol = 1, startRow = cur_row, colNames = FALSE)
    mergeCells(wb, sheet_name, cols = 1:nc, rows = cur_row)
    addStyle(wb, sheet_name,
             createStyle(fontName = "Times New Roman", fontSize = 14,
                         textDecoration = "bold", halign = opt$title_align, valign = "center"),
             rows = cur_row, cols = 1:nc)
    setRowHeights(wb, sheet_name, rows = cur_row, heights = 24)
    cur_row <- cur_row + 1
  }

  # --- 注释 ---
  if (nzchar(note)) {
    writeData(wb, sheet_name, data.frame(x = note),
              startCol = 1, startRow = cur_row, colNames = FALSE)
    mergeCells(wb, sheet_name, cols = 1:nc, rows = cur_row)
    addStyle(wb, sheet_name,
             createStyle(fontName = "Times New Roman", fontSize = 11,
                         fgFill = "#FFFF99", wrapText = TRUE,
                         halign = "left", valign = "top"),
             rows = cur_row, cols = 1:nc)
    note_lines <- max(1, ceiling(nchar(note) / 80))
    setRowHeights(wb, sheet_name, rows = cur_row, heights = note_lines * 14 + 6)
    cur_row <- cur_row + 1
  }

  # --- 表头 + 数据 ---
  header_row <- cur_row
  writeData(wb, sheet_name, sub, startCol = 1, startRow = header_row, colNames = TRUE)
  for (j in seq_len(nc)) {
    addStyle(wb, sheet_name,
             createStyle(fontName = "Times New Roman", textDecoration = "bold",
                         fgFill = palette[[group_of(colnames(sub)[j])]],
                         halign = "center", valign = "center"),
             rows = header_row, cols = j)
  }
  setRowHeights(wb, sheet_name, rows = header_row, heights = 20)

  # --- 数字格式（数据行）：百分比 / 千分位 ---
  if (nrow(sub) > 0) {
    pct_cols   <- which(grepl("Ratio|Rate|Q30|percentage", colnames(sub), ignore.case = TRUE))
    comma_cols <- which(grepl("Reads|Bases", colnames(sub), ignore.case = TRUE))
    data_rows  <- (header_row + 1):(header_row + nrow(sub))
    if (length(pct_cols) > 0)   addStyle(wb, sheet_name, createStyle(numFmt = "0.00%"),
                                         rows = data_rows, cols = pct_cols, gridExpand = TRUE)
    if (length(comma_cols) > 0) addStyle(wb, sheet_name, createStyle(numFmt = "#,##0"),
                                         rows = data_rows, cols = comma_cols, gridExpand = TRUE)
  }

  # --- 列宽（按表头与内容最大长度，上限 40） ---
  for (j in seq_len(nc)) {
    w <- max(nchar(colnames(sub)[j]),
             max(nchar(as.character(sub[[j]])), na.rm = TRUE))
    setColWidths(wb, sheet_name, cols = j, widths = min(w + 2, 40))
  }

  # --- 冻结表头 ---
  freezePane(wb, sheet_name, firstActiveRow = header_row + 1, firstActiveCol = 1)
}

sanitize_sheet <- function(s) {
  s <- as.character(s)
  s <- gsub("[\\[\\]:*?/\\\\]", "_", s)   # Excel 非法 sheet 名字符
  substr(s, 1, 31)                        # Excel sheet 名上限 31 字符
}

sort_df <- function(df) {
  if (nzchar(opt$sort_by) && opt$sort_by %in% colnames(df)) {
    df <- df[order(df[[opt$sort_by]]), , drop = FALSE]
  }
  df
}

# ---------- 差异表多文件合并模式 ----------
if (nzchar(opt$diff_dir)) {
  suffix_re <- paste0(gsub("\\.", "\\\\.", opt$diff_suffix), "$")
  files <- sort(list.files(opt$diff_dir, pattern = suffix_re, full.names = TRUE))
  if (length(files) == 0) stop("diff_dir 下无匹配 *", opt$diff_suffix, " 的文件: ", opt$diff_dir)
  gcol <- opt$group_by
  if (!nzchar(gcol)) stop("差异表模式需 --group_by（如 Regulation）")
  for (f in files) {
    df <- read_df(f)
    if (!(gcol %in% colnames(df))) stop("文件缺少 ", gcol, " 列: ", f)
    df <- sort_df(df)
    filebase <- sub(suffix_re, "", basename(f))
    if (opt$sheet_strip > 0) {
      parts <- strsplit(filebase, "_", fixed = TRUE)[[1]]
      if (length(parts) > opt$sheet_strip) {
        filebase <- paste(parts[1:(length(parts) - opt$sheet_strip)], collapse = "_")
      }
    }
    gvals <- unique(df[[gcol]])
    gvals <- gvals[!is.na(gvals)]   # 跳过 NA（零计数基因 FoldChange=NA → Regulation=NA）
    for (g in gvals) {
      keep <- !is.na(df[[gcol]]) & df[[gcol]] == g
      sub <- df[keep, , drop = FALSE]
      sub[[gcol]] <- NULL
      write_sheet(wb, sanitize_sheet(paste0(g, ".", filebase)), sub,
                  paste0(opt$title, " - ", g), opt$note)
    }
  }
  saveWorkbook(wb, opt$output, overwrite = TRUE)
  cat("已生成 xlsx: ", opt$output, " (", length(files), " 个文件合并)\n", sep = "")
} else {
  df <- read_df(opt$input)
  df <- sort_df(df)
  if (nzchar(opt$group_by) && opt$group_by %in% colnames(df)) {
    gvals <- unique(df[[opt$group_by]])
    gvals <- gvals[!is.na(gvals)]   # 跳过 NA 分组（避免 NA. sheet + 整行 NA 空行）
    for (g in gvals) {
      keep <- !is.na(df[[opt$group_by]]) & df[[opt$group_by]] == g
      sub <- df[keep, , drop = FALSE]
      sub[[opt$group_by]] <- NULL
      gtitle <- if (nzchar(opt$title)) paste0(opt$title, " - ", g) else as.character(g)
      write_sheet(wb, sanitize_sheet(g), sub, gtitle, opt$note)
    }
  } else {
    write_sheet(wb, "Sheet1", df, opt$title, opt$note)
  }
  saveWorkbook(wb, opt$output, overwrite = TRUE)
  cat("已生成 xlsx: ", opt$output, " (", nrow(df), " 行 x ", ncol(df), " 列)\n", sep = "")
}
