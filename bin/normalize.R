#!/usr/bin/env Rscript
# =============================================================================
# normalize.R — Normalize a count matrix to CPM / FPKM / Spike so samples are
# comparable. Each normalization method is a function; which methods to run is
# selected via --methods (comma-separated), default "cpm,fpkm" (both).
#
# Input matrix (tab-separated, header):  gene_id  length  <sample1>  <sample2> ...
# Output matrix:                         gene_id  <s>.<Suffix> ...
#   Suffix per method: Fpkm / Cpm / Spike.
#
#   CPM  / RPM  = count * 1e6 / library_size
#   FPKM / RPKM = count * 1e9 / (length * library_size)
#   Spike       = count * 1e6 / spike_count
# =============================================================================
suppressWarnings(suppressMessages({
    library(argparser)
}))

# ── Parse args ──
argv <- arg_parser("Normalize a count matrix (CPM/FPKM/RPKM, multi-method)")
argv <- add_argument(argv, "--input",   help = "Count matrix (gene_id, length, samples)")
argv <- add_argument(argv, "--methods", help = "Comma-separated: cpm | fpkm | rpkm", default = "cpm,fpkm")
argv <- add_argument(argv, "--total_mapped", help = "Comma-separated *.main_mapped.txt (sample -> pure-main read1 mapped count, spike removed); CPM/FPKM/RPKM denominator")
argv <- add_argument(argv, "--spike_factors", help = "Optional spikein_scale_factors.tsv (sample spike_count factor size_factor); append .Spike column = count x factor (1e6/spike_count), independent of lib_size")
argv <- add_argument(argv, "--output",  help = "Output normalized matrix (tsv)")
argv <- parse_args(argv)

# ── Normalization methods (each takes count matrix + gene length + per-sample denominator) ──
norm_fun <- list(
    cpm  = function(mat, len, lib) sweep(mat, 2, lib, "/") * 1e6,
    fpkm = function(mat, len, lib) sweep(mat, 2, lib, "/") * 1e9 / len,
    rpkm = function(mat, len, lib) sweep(mat, 2, lib, "/") * 1e9 / len   # FPKM alias
)
suffix <- c(fpkm = "Fpkm", cpm = "Cpm", rpkm = "Rpkm")   # method -> column suffix

main <- function(argv) {
    methods <- trimws(unlist(strsplit(argv$methods, ",", fixed = TRUE)))
    methods <- tolower(methods[nzchar(methods)])
    if (length(methods) == 0)
        stop("[normalize] --methods must name at least one method (cpm/fpkm/rpkm)")
    unknown <- setdiff(methods, names(norm_fun))
    if (length(unknown) > 0)
        stop("[normalize] unknown method(s): ", paste(unknown, collapse = ", "),
             " (allowed: cpm, fpkm, rpkm)")

    dt <- read.delim(argv$input, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
                     check.names = FALSE)
    sample_cols <- setdiff(names(dt), c("gene_id", "length"))
    if (length(sample_cols) == 0)
        stop("[normalize] no sample columns found in ", argv$input)

    mat <- as.matrix(dt[, sample_cols, drop = FALSE])
    storage.mode(mat) <- "double"
    len <- as.numeric(dt$length)

    # 归一化分母：--total_mapped 传入每样本 main_mapped（纯主 read1 mapped，已由 EXTRACT_R1 剔 spike）。
    #   不再回退 colSums（基因区计数与 total_mapped 口径不同，混用会导致跨样本 CPM 不可比），
    #   缺失/非正直接报错。
    if (is.null(argv$total_mapped) || is.na(argv$total_mapped) || !nzchar(argv$total_mapped))
        stop("[normalize] --total_mapped 未提供（分母不再回退 colSums）")
    files <- trimws(unlist(strsplit(argv$total_mapped, ",", fixed = TRUE)))
    files <- files[nzchar(files)]
    for (f in files) if (!file.exists(f)) stop("[normalize] total_mapped 文件不存在: ", f)
    lib_size <- setNames(numeric(length(files)), sub("\\.main_mapped\\.txt$", "", basename(files)))
    for (i in seq_along(files)) lib_size[i] <- as.numeric(readLines(files[i], warn = FALSE)[1])
    if (any(!is.finite(lib_size)) || any(lib_size <= 0))
        stop("[normalize] total_mapped 非正或 NA: ",
             paste(names(lib_size)[!is.finite(lib_size) | lib_size <= 0], collapse = ", "))
    miss <- setdiff(sample_cols, names(lib_size))
    if (length(miss) > 0)
        stop("[normalize] 缺 total_mapped 的样本: ", paste(miss, collapse = ", "))
    lib_size <- lib_size[sample_cols]   # 按矩阵列序重排（sweep 按位置对齐，files 收集顺序可能与列序不一致）

    # spike-in：分母已由 EXTRACT_R1 直接产出 main_mapped（纯主，已剔 spike），此处不再减 spike_count。
    #   sf 仅用于下方追加 .Spike 列（count × 1e6/spike_count）。
    sf <- NULL
    if (!is.null(argv$spike_factors) && !is.na(argv$spike_factors) && nzchar(argv$spike_factors)) {
        if (!file.exists(argv$spike_factors))
            stop("[normalize] --spike_factors not found: ", argv$spike_factors)
        sf <- read.delim(argv$spike_factors, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    }
    if (any(lib_size <= 0))
        stop("[normalize] 分母（main_mapped）非正: ",
             paste(names(lib_size)[lib_size <= 0], collapse = ", "))

    # Fixed column order Fpkm -> Cpm -> Rpkm, regardless of --methods order.
    out_order <- c("fpkm", "cpm", "rpkm")
    out_order <- out_order[out_order %in% methods]

    res <- data.frame(gene_id = dt$gene_id, stringsAsFactors = FALSE, check.names = FALSE)
    for (m in out_order) {
        v <- norm_fun[[m]](mat, len, lib_size)
        colnames(v) <- paste0(sample_cols, ".", suffix[[m]])
        res <- cbind(res, v)
    }

    # spike 列：count × factor（factor = 1e6/spike_count，来自 spikein_scale.R），独立追加，不碰 lib_size。
    if (!is.null(sf)) {
        fac <- setNames(as.numeric(sf$factor), sf$sample)
        for (s in sample_cols) {
            if (s %in% names(fac) && !is.na(fac[[s]]) && is.finite(fac[[s]])) {
                res[[paste0(s, ".Spike")]] <- mat[, s] * fac[[s]]
            } else {
                message("[normalize] WARNING: no valid spike factor for ", s, " — .Spike omitted")
            }
        }
        message("[normalize] appended .Spike columns (count x spike factor)")
    }

    write.table(res, file = argv$output, sep = "\t", quote = FALSE, row.names = FALSE)
    message("[normalize] methods=", paste(out_order, collapse = ","), " -> ", argv$output)
}

main(argv)
