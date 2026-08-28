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
# =============================================================================
suppressWarnings(suppressMessages({
    library(argparser)
}))

# ── Parse args ──
argv <- arg_parser("Normalize a count matrix (CPM/FPKM/RPKM, multi-method)")
argv <- add_argument(argv, "--input",   help = "Count matrix (gene_id, length, samples)")
argv <- add_argument(argv, "--methods", help = "Comma-separated: cpm | fpkm | rpkm", default = "cpm,fpkm")
argv <- add_argument(argv, "--output",  help = "Output normalized matrix (tsv)")
argv <- parse_args(argv)

# ── Normalization methods (each takes count matrix + gene length) ──
norm_fun <- list(
    cpm  = function(mat, len) sweep(mat, 2, colSums(mat, na.rm = TRUE), "/") * 1e6,
    fpkm = function(mat, len) sweep(mat, 2, colSums(mat, na.rm = TRUE), "/") * 1e9 / len,
    rpkm = function(mat, len) sweep(mat, 2, colSums(mat, na.rm = TRUE), "/") * 1e9 / len   # FPKM alias
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

    lib_size <- colSums(mat, na.rm = TRUE)
    if (any(lib_size == 0))
        message("[normalize] warning: zero library size for sample(s): ",
                paste(names(lib_size)[lib_size == 0], collapse = ", "))

    # Fixed column order Fpkm -> Cpm -> Rpkm, regardless of --methods order.
    out_order <- c("fpkm", "cpm", "rpkm")
    out_order <- out_order[out_order %in% methods]

    res <- data.frame(gene_id = dt$gene_id, stringsAsFactors = FALSE, check.names = FALSE)
    for (m in out_order) {
        v <- norm_fun[[m]](mat, len)
        colnames(v) <- paste0(sample_cols, ".", suffix[[m]])
        res <- cbind(res, v)
    }

    write.table(res, file = argv$output, sep = "\t", quote = FALSE, row.names = FALSE)
    message("[normalize] methods=", paste(out_order, collapse = ","), " -> ", argv$output)
}

main(argv)
