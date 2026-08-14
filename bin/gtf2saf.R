#!/usr/bin/env Rscript
# =============================================================================
# gtf2saf.R — Generate TSS and gene body SAF annotations from GTF
#
# Reads GTF once, outputs two SAF files:
#   1. TSS window:      TSS-tss_upstream to TSS+tss_downstream
#   2. Gene body:       TSS+genebody_offset    to TES
#
# SAF Format (featureCounts requires these exact column names, case-sensitive)
# GeneID    Chr Start   End Strand
# GeneA chr1    100 400 +
# GeneB chr2    1500    2000    -
#
# Plus-strand genes:  TSS = start, TES = end
# Minus-strand genes: TSS = end,   TES = start
# =============================================================================
suppressWarnings(suppressMessages({
    library(argparser)
    library(rtracklayer)
}))

# ── Parse args ──
argv <- arg_parser("Generate TSS window and gene body SAF annotations from GTF")
argv <- add_argument(argv, "--gtf",              help = "Path to reference GTF file")
argv <- add_argument(argv, "--tss_upstream",     help = "TSS upstream window (bp)",    default = 50,  type = "integer")
argv <- add_argument(argv, "--tss_downstream",   help = "TSS downstream window (bp)",  default = 300, type = "integer")
argv <- add_argument(argv, "--genebody_offset",  help = "Gene body start offset (bp)", default = 301, type = "integer")
argv <- add_argument(argv, "--outdir",           help = "Output dir", default = "./")
argv <- parse_args(argv)

# ── Read GTF ──
read_gtf <- function(gtf_file) {
    message("[gtf2saf] Reading GTF: ", gtf_file)
    gr <- rtracklayer::import(gtf_file, format = "gtf")   # GRanges; attributes parsed into mcols
    gr <- gr[gr$type == "gene"]                           # keep gene-level entries only
    if (is.null(mcols(gr)$gene_id)) stop("GTF missing 'gene_id' attribute")
    message("[gtf2saf] ", length(gr), " gene entries")
    return(gr)
}

# ── Validate inputs ──
check_inputs <- function(argv) {
    if (!file.exists(argv$gtf)) stop("GTF file not found: ", argv$gtf)
    if (argv$genebody_offset <= argv$tss_downstream)
        stop("genebody_offset (", argv$genebody_offset, ") must be > tss_downstream (",
             argv$tss_downstream, ") so the gene body does not overlap the TSS window")
}

# ── Build SAF data.frame from per-strand coordinate vectors ──
build_saf <- function(gene_ids, chr, strand, plus_start, plus_end, minus_start, minus_end) {
    saf_start <- saf_end <- integer(length(gene_ids))
    plus  <- strand == "+"
    minus <- strand == "-"

    saf_start[plus]  <- plus_start[plus]
    saf_end[plus]    <- plus_end[plus]
    saf_start[minus] <- minus_start[minus]
    saf_end[minus]   <- minus_end[minus]

    data.frame(
        GeneID = gene_ids, Chr = chr,
        Start  = pmax(saf_start, 1), End = saf_end,
        Strand = strand, stringsAsFactors = FALSE
    )
}

# ── Main ──
main <- function(argv) {
    check_inputs(argv)

    tss_upstream     <- as.integer(argv$tss_upstream)
    tss_downstream   <- as.integer(argv$tss_downstream)
    genebody_offset  <- as.integer(argv$genebody_offset)

    gr <- read_gtf(argv$gtf)

    gene_ids <- mcols(gr)$gene_id
    chr    <- as.character(seqnames(gr))
    start  <- start(gr)
    end    <- end(gr)
    strand <- as.character(strand(gr))

    # ── Report genes too short for a valid gene body ──
    gene_len <- end - start + 1
    n_short  <- sum(gene_len < genebody_offset)
    message("[gtf2saf] genes shorter than genebody_offset (", genebody_offset,
            " bp): ", n_short, "/", length(gene_len))

    # ── Build TSS SAF ──
    tss_saf <- build_saf(
        gene_ids, chr, strand,
        plus_start  = start - tss_upstream,
        plus_end    = pmin(start + tss_downstream, end),   # clamp to gene end (TES)
        minus_start = pmax(end   - tss_downstream, start), # clamp to gene start (TES)
        minus_end   = end   + tss_upstream
    )

    # ── Build gene body SAF ──
    gb_saf <- build_saf(
        gene_ids, chr, strand,
        plus_start  = start + genebody_offset, plus_end  = end,
        minus_start = start,                    minus_end = end - genebody_offset
    )

    # ── Filter invalid ──
    tss_ok <- !is.na(tss_saf$GeneID) & tss_saf$GeneID != "" & tss_saf$Start < tss_saf$End
    gb_ok  <- !is.na(gb_saf$GeneID)  & gb_saf$GeneID  != "" & gb_saf$Start  < gb_saf$End

    tss_saf <- tss_saf[tss_ok, ]
    gb_saf  <- gb_saf[gb_ok, ]

    # ── Deduplicate ──
    tss_saf <- tss_saf[!duplicated(tss_saf$GeneID), ]
    gb_saf  <- gb_saf[!duplicated(gb_saf$GeneID), ]

    # ── Write outputs ──
    write.table(tss_saf, file = file.path(argv$outdir, "tss.saf"), sep = "\t", quote = FALSE, row.names = FALSE)
    message("[gtf2saf] TSS SAF: ", nrow(tss_saf), " regions → ", file.path(argv$outdir, "tss.saf"))

    write.table(gb_saf, file = file.path(argv$outdir, "genebody.saf"), sep = "\t", quote = FALSE, row.names = FALSE)
    message("[gtf2saf] Gene body SAF: ", nrow(gb_saf), " regions → ", file.path(argv$outdir, "genebody.saf"))

    message("[gtf2saf] Done.")
}

main(argv)
