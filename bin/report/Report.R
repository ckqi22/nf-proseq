#!/usr/bin/env Rscript
# =============================================================================
# report.R
# =============================================================================
suppressPackageStartupMessages({
  library(argparse)
  library(openxlsx)
  library(yaml)
  library(dplyr)
  library(reticulate)
})

########################## command line arguments ##########################
parser <- ArgumentParser(description = "Package PRO-seq results into a client deliverable folder")
parser$add_argument('-o', '--output_dir', default = '.', help = 'output directory')
parser$add_argument('--config', default = 'params.yml', help = 'params.yml config file')
parser$add_argument('--samplesheet', default = 'samplesheet.csv', help = 'samplesheet csv file')
parser$add_argument('--diff_dir', help = 'DESeq2 output directory (deprecated; xlsx now from --xlsx_dir)')
parser$add_argument('--enrich_dir', help = 'GO/KEGG/GSEA enrichment results directory')
parser$add_argument('--metagene_dir', help = 'metagene results directory')
parser$add_argument('--pausing_dir', help = 'pausing index results directory')
parser$add_argument('--pol2_signal', help = 'pol2_signal_table.tsv (Pol II active site)')
parser$add_argument('--pol2_signal_xlsx', help = 'PROSeq_pol2_signal.xlsx (Pol II active site, with note)')
parser$add_argument('--plot_dir', help = 'scatter/volcano/heatmap figure directory')
parser$add_argument('--statistics', help = 'read statistics file (deprecated; xlsx now from --xlsx_dir)')
parser$add_argument('--xlsx_dir', help = 'txt2xlsx output directory (TXT2XLSX 模块产物)')
parser$add_argument('--resolved_config', help = 'config.txt from parse_config.py (for build/gtf resolution)')
parser$add_argument('--fastqc_images', help = 'extract_fastqc_images.sh output dir (raw/ + trimmed/)')
parser$add_argument('--base_quality_plot', help = 'base_quality_plot.R output dir (*.tiff/*.pdf: read1/2_before/after_filtering)')
parser$add_argument('--group_heatmaps', help = 'group-level metagene heatmap dir (from metagene_group)')
parser$add_argument('--resources_dir', default = '/workplace/pipeline/resources', help = 'resources directory')

args <- parser$parse_args()

########################## normalize paths ##########################
output_dir   <- normalizePath(gsub('/$', '', args$output_dir), mustWork = FALSE)
config_file  <- normalizePath(args$config)
samplesheet  <- normalizePath(args$samplesheet)
resources_dir <- normalizePath(gsub('/$', '', args$resources_dir))

# 可选路径统一走 normalize_optional：缺失即 NULL，绝不因路径不存在而中断
normalize_optional <- function(p) {
  if (is.null(p) || is.na(p) || nchar(p) == 0) return(NULL)
  normalizePath(p, mustWork = FALSE)
}
diff_dir        <- normalize_optional(args$diff_dir)
enrich_dir      <- normalize_optional(args$enrich_dir)
metagene_dir    <- normalize_optional(args$metagene_dir)
pausing_dir     <- normalize_optional(args$pausing_dir)
plot_dir        <- normalize_optional(args$plot_dir)
pol2_signal     <- normalize_optional(args$pol2_signal)
pol2_signal_xlsx <- normalize_optional(args$pol2_signal_xlsx)
statistics      <- normalize_optional(args$statistics)
xlsx_dir        <- normalize_optional(args$xlsx_dir)
resolved_config <- normalize_optional(args$resolved_config)
fastqc_images   <- normalize_optional(args$fastqc_images)
base_quality_plot <- normalize_optional(args$base_quality_plot)
group_heatmaps  <- normalize_optional(args$group_heatmaps)

########################## read config (params.yml) ##########################
custom_handlers <- list(
  'bool#yes' = function(x) { if (x == 'y' || x == 'Y') x else if (x == 'n' || x == 'N') x else TRUE },
  'bool#no'  = function(x) { if (x == 'n' || x == 'N') x else FALSE }
)
config <- yaml.load(readLines(config_file), handlers = custom_handlers)

get_field <- function(config, key, default = '') {
  v <- config[[key]]
  if (is.null(v) || length(v) == 0 || (length(v) == 1 && is.na(v))) return(default)
  as.character(v[1])
}

name          <- get_field(config, 'name', 'PRO-seq')
project_no    <- get_field(config, 'projectNo.')
institute     <- get_field(config, 'institute')
species       <- get_field(config, 'species')
build         <- get_field(config, 'build')
instrument    <- get_field(config, 'instrument_model')
sample_type   <- get_field(config, 'sample_type')
sample_number <- get_field(config, 'sample_number')
strandedness  <- get_field(config, 'strandedness', 'reverse')
gtf_file      <- get_field(config, 'gtf')
species_title <- gsub(' ', '_', tools::toTitleCase(tolower(species)))
current_date  <- paste(unlist(strsplit(as.character(Sys.Date()), '-')), collapse = '')

########################## build/gtf 解析（复用 parse_config.py 的 config.txt） ##########################
# params.yml 的 build/gtf 常为空；此时从 --resolved_config（parse_config.py 产出的
# config.txt，含 `build:`/`gtf:` 行）解析，使 IGV 章节能拷到 gtf。
resolved_build <- NULL
resolved_gtf   <- NULL
if (!is.null(resolved_config) && file.exists(resolved_config)) {
  lines <- readLines(resolved_config, warn = FALSE)
  for (ln in lines) {
    kv <- strsplit(ln, ':', fixed = TRUE)[[1]]
    if (length(kv) < 2) next
    k <- trimws(kv[1])
    v <- trimws(paste(kv[-1], collapse = ':'))
    if (k == 'build') resolved_build <- v
    else if (k == 'gtf') resolved_gtf <- v
  }
}
if (is.null(build) || build == '' || build == 'null' || is.na(build)) build <- resolved_build
if (is.null(gtf_file) || gtf_file == '' || is.na(gtf_file)) gtf_file <- resolved_gtf
if (is.null(build) || is.na(build)) build <- ''

output_dir <- file.path(output_dir,
                        paste0(name, '_', project_no, '_', species_title,
                               '_PRO-seq_Sequencing_Report_', current_date))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

########################## helpers ##########################
# 安全 cp：glob 无匹配时打印 [skip] 并跳过（不再静默吞错）
safe_cp <- function(pattern, dest_dir, recursive = FALSE) {
  files <- Sys.glob(pattern)
  if (length(files) == 0) {
    cat('[skip] 无匹配文件:', pattern, '\n')
    return(invisible(NULL))
  }
  flag <- if (recursive) ' -r' else ''
  for (f in files) {
    if (!file.exists(f)) next
    system(paste0('cp', flag, ' ', shQuote(f), ' ', shQuote(dest_dir)))
  }
  invisible(files)
}

# 单文件拷贝到指定路径（可改名）；源缺失即 [skip]
safe_cp_to <- function(src, dest_path) {
  if (is.null(src) || !file.exists(src)) {
    cat('[skip] 无匹配文件:', if (is.null(src)) '<NULL>' else src, '\n')
    return(invisible(FALSE))
  }
  file.copy(src, dest_path, overwrite = TRUE)
  invisible(TRUE)
}

# 富集结果目录统一命名（enrichment.R 产出的 go/pathway/gsea 子目录 → 报告规范名）
#   go.up.Gene.A_vs_B_GeneBody_paired        -> up.go.gene.A_vs_B
#   pathway.down.Gene.A_vs_B_GeneBody_unpaired -> down.pathway.gene.A_vs_B
#   gsea.Gene.A_vs_B_GeneBody_paired.pathway   -> gsea.Gene.A_vs_B.pathway（保留区分后缀）
rename_enrich <- function(name) {
  if (grepl('^go\\.(up|down)\\.Gene\\.', name))
    return(sub('^go\\.(up|down)\\.Gene\\.(.+)_GeneBody_(paired|unpaired)$', '\\1.go.gene.\\2', name))
  if (grepl('^pathway\\.(up|down)\\.Gene\\.', name))
    return(sub('^pathway\\.(up|down)\\.Gene\\.(.+)_GeneBody_(paired|unpaired)$', '\\1.pathway.gene.\\2', name))
  if (grepl('^gsea\\.Gene\\.', name))
    return(sub('^gsea\\.Gene\\.(.+)_GeneBody_(paired|unpaired)\\.(pathway|go)$', 'gsea.Gene.\\1.\\3', name))
  name
}

# glob 匹配子目录并按 rename_fun 改名后 cp -r 到 dest_dir
safe_cp_rename <- function(pattern, dest_dir, rename_fun) {
  files <- Sys.glob(pattern)
  if (length(files) == 0) {
    cat('[skip] 无匹配目录:', pattern, '\n')
    return(invisible(NULL))
  }
  for (f in files) {
    if (!dir.exists(f)) next
    system(paste0('cp -r ', shQuote(f), ' ', shQuote(file.path(dest_dir, rename_fun(basename(f))))))
  }
  invisible(files)
}

# 通用样式
title_style  <- createStyle(fontSize = 16, textDecoration = 'Bold')
common_style <- createStyle(fontSize = 11, textDecoration = 'Bold', fgFill = '#AECDD7',
                            halign = 'center', valign = 'center')
centre_style <- createStyle(halign = 'center', valign = 'center')

########################## 01. Project_Info ##########################
sample_information <- function(samplesheet, config, build, output_dir) {
  ss <- read.table(samplesheet, header = TRUE, sep = ',', check.names = FALSE,
                   fill = TRUE, stringsAsFactors = FALSE, comment.char = '')
  samples <- ss$sample
  groups  <- if ('group' %in% colnames(ss)) ss$group else rep('unknown', nrow(ss))
  groups[is.na(groups) | groups == ''] <- 'unknown'

  sample_df <- data.frame(
    `Sample ID`      = seq_along(samples),
    `Sample Name`    = samples,
    `Group Name`     = groups,
    `Quality Status` = rep('OK', length(samples)),
    check.names      = FALSE
  )

  wb <- createWorkbook()
  addWorksheet(wb, 'Sheet1')
  setColWidths(wb, 'Sheet1', cols = 1:4, widths = 20)
  modifyBaseFont(wb, fontSize = 11, fontColour = 'black', fontName = 'Times New Roman')

  writeData(wb, 'Sheet1', 'Table 1. Sample Information', startRow = 1, startCol = 1)
  info_rows <- list(
    c('Species', tools::toTitleCase(tolower(species))),
    c('Sample type', sample_type),
    c('Sample number', sample_number),
    c('Instrument model', instrument),
    c('Genome build', build),
    c('Library orientation (R1 vs nascent RNA)', strandedness)
  )
  for (i in seq_along(info_rows)) {
    writeData(wb, 'Sheet1', info_rows[[i]][1], startRow = 1 + i, startCol = 1)
    writeData(wb, 'Sheet1', info_rows[[i]][2], startRow = 1 + i, startCol = 2)
  }

  data_row <- 8
  writeData(wb, 'Sheet1', sample_df, startRow = data_row, startCol = 1, rowNames = FALSE)

  addStyle(wb, 'Sheet1', title_style, rows = 1, cols = 1)
  addStyle(wb, 'Sheet1', common_style, rows = 2:(data_row - 1), cols = 1)
  addStyle(wb, 'Sheet1', common_style, rows = data_row, cols = 1:4)
  addStyle(wb, 'Sheet1', centre_style, rows = (data_row + 1):(data_row + nrow(sample_df)),
           cols = 1:4, gridExpand = TRUE)

  saveWorkbook(wb, file.path(output_dir, 'Project_Info.xlsx'), overwrite = TRUE)
}

########################## run sections ##########################
### mkdir
dirs <- c('01.Project_Info', '02.Lab_QC', '03.Data_QC/Mean_Base_Quality',
          '04.nascent_mRNA_Profiling', '05.Diff_nascent_mRNAs',
          '06.nascent_mRNA_GO', '07.nascent_mRNA_Pathway', '08.nascent_mRNA_GSEA',
          '09.Cluster', '10.Scatter', '11.Volcano',
          '12.Meta_TSS', '13.Pausing_analysis',
          '14.Pol_II_active_site/14.1.PROSeq_profiling', '14.Pol_II_active_site/14.2.Heatmap',
          '15.IGV', '16.GEO')
for (d in dirs) dir.create(file.path(output_dir, d), recursive = TRUE, showWarnings = FALSE)

# 01. 样本信息（唯一本脚本生成的 xlsx）
sample_information(samplesheet, config, build, file.path(output_dir, '01.Project_Info'))

# 03. Data_QC：Read_statistics.xlsx（改名）+ Mean_Base_Quality（过滤前/后碱基质量分布图）
if (!is.null(xlsx_dir)) {
  safe_cp_to(file.path(xlsx_dir, 'read_statistics.xlsx'),
             file.path(output_dir, '03.Data_QC', 'Read_statistics.xlsx'))
} else {
  cat('[skip] 未指定 --xlsx_dir，跳过 Read_statistics\n')
}
if (!is.null(base_quality_plot) && dir.exists(base_quality_plot)) {
  safe_cp(file.path(base_quality_plot, '*.tiff'), file.path(output_dir, '03.Data_QC', 'Mean_Base_Quality'))
  safe_cp(file.path(base_quality_plot, '*.pdf'),  file.path(output_dir, '03.Data_QC', 'Mean_Base_Quality'))
} else {
  cat('[skip] 未指定 --base_quality_plot 或目录不存在\n')
}

# 04. 表达谱（txt2xlsx 产出，改名为 nascent_mRNA_Profiling.xlsx）
if (!is.null(xlsx_dir)) {
  safe_cp_to(file.path(xlsx_dir, 'expression_profiling.xlsx'),
             file.path(output_dir, '04.nascent_mRNA_Profiling', 'nascent_mRNA_Profiling.xlsx'))
}

# 05. 差异表（txt2xlsx 产出）
if (!is.null(xlsx_dir)) {
  safe_cp(file.path(xlsx_dir, 'All_Comparisons_genes.xlsx'),          file.path(output_dir, '05.Diff_nascent_mRNAs'))
  safe_cp(file.path(xlsx_dir, 'Differentially_Expressed_genes.xlsx'), file.path(output_dir, '05.Diff_nascent_mRNAs'))
}

# 06 / 07 / 08 富集（go.* / pathway.* / gsea.* 子目录统一改名）+ 静态 Summary
if (!is.null(enrich_dir) && dir.exists(enrich_dir)) {
  safe_cp_rename(file.path(enrich_dir, 'go.*'),     file.path(output_dir, '06.nascent_mRNA_GO'),       rename_enrich)
  safe_cp_rename(file.path(enrich_dir, 'pathway.*'), file.path(output_dir, '07.nascent_mRNA_Pathway'), rename_enrich)
  safe_cp_rename(file.path(enrich_dir, 'gsea.*'),    file.path(output_dir, '08.nascent_mRNA_GSEA'),    rename_enrich)
}
safe_cp(file.path(resources_dir, 'report', 'Summary_GO.pdf'),     file.path(output_dir, '06.nascent_mRNA_GO'))
safe_cp(file.path(resources_dir, 'report', 'Summary_PATHWAY.pdf'), file.path(output_dir, '07.nascent_mRNA_Pathway'))

# 09 / 10 / 11 图 + 说明
if (!is.null(plot_dir) && dir.exists(plot_dir)) {
  safe_cp(file.path(plot_dir, 'heatmap.*'), file.path(output_dir, '09.Cluster'))
  safe_cp(file.path(plot_dir, 'scatter.*'), file.path(output_dir, '10.Scatter'))
  safe_cp(file.path(plot_dir, 'volcano.*'), file.path(output_dir, '11.Volcano'))
}
safe_cp(file.path(resources_dir, 'report', '热图说明.txt'), file.path(output_dir, '09.Cluster'))
safe_cp(file.path(resources_dir, 'report', '散点图说明.txt'), file.path(output_dir, '10.Scatter'))
safe_cp(file.path(resources_dir, 'report', '火山图说明.txt'), file.path(output_dir, '11.Volcano'))

# 13. Pausing_analysis（仅 PROSeq_pausing.xlsx）
if (!is.null(pausing_dir) && dir.exists(pausing_dir)) {
  safe_cp(file.path(pausing_dir, 'PROSeq_pausing.xlsx'), file.path(output_dir, '13.Pausing_analysis'))
} else {
  cat('[skip] 未指定 --pausing_dir 或目录不存在\n')
}

# 12. Meta_TSS（仅 *metagene_profile.pdf + *.tiff；.tiff 由 PLOTPROFILE 管线内转出）
if (!is.null(metagene_dir) && dir.exists(metagene_dir)) {
  safe_cp(file.path(metagene_dir, '*_metagene_profile_*.pdf'),  file.path(output_dir, '12.Meta_TSS'))
  safe_cp(file.path(metagene_dir, '*_metagene_profile_*.tiff'), file.path(output_dir, '12.Meta_TSS'))
} else {
  cat('[skip] 未指定 --metagene_dir 或目录不存在\n')
}

# 14. Pol_II_active_site（14.1.PROSeq_profiling 仅交付 xlsx + 14.2.Heatmap 组级热图）
if (!is.null(pol2_signal_xlsx) && file.exists(pol2_signal_xlsx)) {
  safe_cp(pol2_signal_xlsx, file.path(output_dir, '14.Pol_II_active_site', '14.1.PROSeq_profiling'))
} else {
  cat('[skip] pol2_signal_xlsx 缺失（14.1 留空占位）\n')
}
if (!is.null(group_heatmaps) && dir.exists(group_heatmaps)) {
  safe_cp(file.path(group_heatmaps, '*_metagene_heatmap_*.pdf'),  file.path(output_dir, '14.Pol_II_active_site', '14.2.Heatmap'))
  safe_cp(file.path(group_heatmaps, '*_metagene_heatmap_*.tiff'), file.path(output_dir, '14.Pol_II_active_site', '14.2.Heatmap'))
} else {
  cat('[skip] 未指定 --group_heatmaps 或目录不存在\n')
}

# 15 IGV
if (!is.null(gtf_file) && nchar(gtf_file) > 0 && file.exists(gtf_file)) {
  system(paste0('cp ', shQuote(gtf_file), ' ', shQuote(file.path(output_dir, '15.IGV'))))
  system(paste0('gzip -f ', shQuote(file.path(output_dir, '15.IGV', basename(gtf_file)))))
} else {
  cat('[skip] gtf 缺失（params.yml 的 gtf 为空且 --resolved_config 未解析到），跳过 IGV gtf\n')
}
safe_cp(file.path(resources_dir, 'report', 'IGV_2.16.2.zip'), file.path(output_dir, '15.IGV'))
safe_cp(file.path(resources_dir, 'report', '云序生物超详细_IGV可视化操作详解.pdf'), file.path(output_dir, '15.IGV'))
safe_cp(file.path(resources_dir, 'report', 'IGV可视化说明.txt'), file.path(output_dir, '15.IGV'))

# 16 GEO
safe_cp(file.path(output_dir, '04.nascent_mRNA_Profiling', 'nascent_mRNA_Profiling.xlsx'), file.path(output_dir, '16.GEO'))
safe_cp(file.path(resources_dir, 'report', 'FileZilla_3.11.0.2_win64-setup.exe'), file.path(output_dir, '16.GEO'))
safe_cp(file.path(resources_dir, 'report', 'FileZilla_3.12.0.2_win32-setup.exe'), file.path(output_dir, '16.GEO'))
safe_cp(file.path(resources_dir, 'report', '云序生物GEO原始数据上传流程说明.pdf'), file.path(output_dir, '16.GEO'))
safe_cp(file.path(resources_dir, 'report', 'seq_template_注释版_.xlsx'), file.path(output_dir, '16.GEO'))

########################## Word 模板盖章 ##########################
word_template <- file.path(resources_dir, 'word', '云序生物微量PRO-seq测序报告.docx')
word_save     <- file.path(output_dir, paste0(name, '_PRO-seq_测序报告.docx'))

stamp_docx <- function() {
  if (!file.exists(word_template)) {
    message('skip docx stamping: template not found: ', word_template)
    return(invisible(FALSE))
  }
  env_path <- '/workplace/hanguojun/mambaforge/envs/snakemake'
  python_path <- file.path(env_path, 'bin/python')
  if (file.exists(python_path)) {
    try(reticulate::use_python(python_path, required = FALSE), silent = TRUE)
  }
  reticulate::source_python('/workplace/pipeline/code/modify_docx.py')
  modify_doc(word_template, word_save, name, institute, project_no)
  invisible(TRUE)
}

stamp_ok <- tryCatch({ stamp_docx(); TRUE }, error = function(e) {
  message('docx stamping failed (non-fatal): ', conditionMessage(e))
  FALSE
})

cat(paste0('PRO-seq report generated: ', output_dir, '\n'))
if (stamp_ok) cat('docx stamped\n') else cat('docx stamping skipped/failed\n')
