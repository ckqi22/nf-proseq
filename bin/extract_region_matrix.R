#!/usr/bin/env Rscript
# =============================================================================
# extract_region_matrix.R — 从 merged counts 文件中提取指定区域的样本矩阵
#
# 输入：多个 *_proseq_counts.txt 文件（每文件含 full_gene / tss / gene_body 三列）
# 输出：合并后的矩阵（gene_id 行 × 样本列）
# =============================================================================
suppressWarnings(suppressMessages({
    library(argparser)
}))

argv <- arg_parser("Extract region-specific columns from proseq merged count files")
argv <- add_argument(argv, "--files",  help = "Comma-separated list of input files")
argv <- add_argument(argv, "--region", help = "Region suffix: _full_gene, _tss, or _gene_body")
argv <- add_argument(argv, "--output", help = "Output merged matrix file")
argv <- parse_args(argv)

files  <- strsplit(argv$files, ",")[[1]]
region <- argv$region

message("[extract_region] Region: ", region, " | ", length(files), " files")

first <- read.delim(files[1], header = TRUE, stringsAsFactors = FALSE)

# Find all columns ending with the region suffix
target_cols <- grep(paste0(region, "$"), names(first), value = TRUE)
if (length(target_cols) == 0) {
    stop("No columns found with suffix: ", region, " in ", files[1])
}

# Build result: gene_id + first file's target column
result <- first[, c("gene_id", target_cols[1]), drop = FALSE]

# Add remaining files
for (i in seq_along(files)[-1]) {
    dt  <- read.delim(files[i], header = TRUE, stringsAsFactors = FALSE)
    tc  <- grep(paste0(region, "$"), names(dt), value = TRUE)
    if (length(tc) == 0) next
    result[[tc[1]]] <- dt[[tc[1]]][match(result$gene_id, dt$gene_id)]
}

result[is.na(result)] <- 0

message("[extract_region] Output: ", nrow(result), " genes x ", ncol(result) - 1, " samples")
write.table(result, file = argv$output, sep = "\t", quote = FALSE, row.names = FALSE)
message("[extract_region] Done: ", argv$output)
