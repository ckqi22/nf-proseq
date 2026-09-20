#!/usr/bin/env Rscript
# =============================================================================
# normalize.R — Normalize a count matrix to CPM / FPKM / RPKM so samples are
# comparable. Each normalization method is a function; which methods to run is
# selected via --methods (comma-separated), default "cpm,fpkm" (both).
#
# Input matrix (tab-separated, header):  gene_id  length  <sample1>  <sample2> ...
# Output matrix:                         gene_id  <s>.<Suffix> ...
#   Suffix per method: Fpkm / Cpm / Rpkm. Length and raw counts are dropped.
#
#   CPM  / RPM  = count / library_size * 1e6
#   FPKM / RPKM = count * 1e9 / (length * library_size)
#   library_size（回退链）：total_mapped > colSums（.Cpm/.Fpkm/.Rpkm 分母恒用 read1 mapped 总数）。
#     --total_mapped 给定时换成每样本 read1 mapped 总数；--spike_factors 给定时额外追加 .Spike 列
#     （= count × factor，factor=1e6/spike_count），不碰 lib_size。
# =============================================================================
suppressWarnings(suppressMessages({
    library(argparser)
}))

# ── Parse args ──
argv <- arg_parser("Normalize a count matrix (CPM/FPKM/RPKM, multi-method)")
argv <- add_argument(argv, "--input",   help = "Count matrix (gene_id, length, samples)")
argv <- add_argument(argv, "--methods", help = "Comma-separated: cpm | fpkm | rpkm", default = "cpm,fpkm")
argv <- add_argument(argv, "--total_mapped", help = "Comma-separated *.total_mapped.txt (sample -> read1 mapped count); use as denominator instead of colSums")
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

    # 归一化分母：total_mapped > colSums（.Cpm/.Fpkm/.Rpkm 恒用 read1 mapped 总数，末级回退 colSums）。
    #   spike 归一化独立成 .Spike 列（见下），不覆盖 lib_size。
    lib_size <- colSums(mat, na.rm = TRUE)
    if (!is.null(argv$total_mapped) && !is.na(argv$total_mapped) && nzchar(argv$total_mapped)) {
        files <- trimws(unlist(strsplit(argv$total_mapped, ",", fixed = TRUE)))
        files <- files[nzchar(files)]
        tm <- setNames(numeric(length(files)), sub("\\.total_mapped\\.txt$", "", basename(files)))
        for (i in seq_along(files)) tm[i] <- as.numeric(readLines(files[i], warn = FALSE)[1])
        for (s in sample_cols)
            if (s %in% names(tm) && !is.na(tm[[s]]) && tm[[s]] > 0) lib_size[[s]] <- tm[[s]]
    }
    if (any(lib_size == 0))
        stop("[normalize] zero library size for sample(s): ",
             paste(names(lib_size)[lib_size == 0], collapse = ", "),
             " — cannot normalize (division by zero)")

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
    if (!is.null(argv$spike_factors) && !is.na(argv$spike_factors) && nzchar(argv$spike_factors)) {
        if (!file.exists(argv$spike_factors))
            stop("[normalize] --spike_factors not found: ", argv$spike_factors)
        sf <- read.delim(argv$spike_factors, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
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
