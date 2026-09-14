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
#   --spike_factors 给定时，library_size 换成每样本 spike_count（spike-in 缩放）。
# =============================================================================
suppressWarnings(suppressMessages({
    library(argparser)
}))

# ── Parse args ──
argv <- arg_parser("Normalize a count matrix (CPM/FPKM/RPKM, multi-method)")
argv <- add_argument(argv, "--input",   help = "Count matrix (gene_id, length, samples)")
argv <- add_argument(argv, "--methods", help = "Comma-separated: cpm | fpkm | rpkm", default = "cpm,fpkm")
argv <- add_argument(argv, "--spike_factors", help = "Optional spikein_scale_factors.tsv (sample spike_count factor size_factor); use spike_count as denominator instead of colSums")
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

    # 归一化分母：默认库大小（colSums）；给了 --spike_factors 时，用每样本 spike_count 替换
    # （spike 缩放 = count/spike_count*1e6，跨样本可比）。spike_count 缺失或 <=0 的样本退回库大小。
    lib_size <- colSums(mat, na.rm = TRUE)
    if (!is.null(argv$spike_factors) && !is.na(argv$spike_factors) && nzchar(argv$spike_factors)) {
        if (!file.exists(argv$spike_factors))
            stop("[normalize] --spike_factors not found: ", argv$spike_factors)
        sf <- read.delim(argv$spike_factors, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
        spk <- setNames(as.numeric(sf$spike_count), sf$sample)
        for (s in sample_cols) {
            if (s %in% names(spk) && !is.na(spk[[s]]) && spk[[s]] > 0) {
                lib_size[[s]] <- spk[[s]]
            } else {
                message("[normalize] WARNING: no valid spike_count for ", s,
                        " — fallback to library size")
            }
        }
        message("[normalize] using spike-in counts as normalization denominator")
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

    write.table(res, file = argv$output, sep = "\t", quote = FALSE, row.names = FALSE)
    message("[normalize] methods=", paste(out_order, collapse = ","), " -> ", argv$output)
}

main(argv)
