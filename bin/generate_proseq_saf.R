#!/usr/bin/env Rscript
# =============================================================================
# generate_proseq_saf.R — Generate TSS and gene body SAF annotations from GTF
#
# Reads GTF once, outputs two SAF files:
#   1. TSS window:      TSS-upstream to TSS+downstream
#   2. Gene body:       TSS+offset    to TES
#
# Plus-strand genes:  TSS = start, TES = end
# Minus-strand genes: TSS = end,   TES = start
# =============================================================================
suppressWarnings(suppressMessages({
    library(argparser)
}))

argv <- arg_parser("Generate TSS window and gene body SAF annotations from GTF")
argv <- add_argument(argv, "--gtf",        help = "Path to reference GTF file")
argv <- add_argument(argv, "--upstream",   help = "TSS upstream window (bp)",    default = 50)
argv <- add_argument(argv, "--downstream", help = "TSS downstream window (bp)",  default = 300)
argv <- add_argument(argv, "--offset",     help = "Gene body start offset (bp)", default = 301)
argv <- add_argument(argv, "--tss_out",    help = "Output TSS SAF file path")
argv <- add_argument(argv, "--gb_out",     help = "Output gene body SAF file path")
argv <- parse_args(argv)

upstream   <- as.integer(argv$upstream)
downstream <- as.integer(argv$downstream)
offset     <- as.integer(argv$offset)

if (!file.exists(argv$gtf)) stop("GTF file not found: ", argv$gtf)

# ── Read GTF ──
message("[generate_saf] Reading GTF: ", argv$gtf)
gtf <- read.delim(argv$gtf, header = FALSE, stringsAsFactors = FALSE, comment.char = "#")
names(gtf) <- c("chr", "source", "feature", "start", "end", "score", "strand", "frame", "attributes")

# Use gene-level entries
gtf <- gtf[gtf$feature == "gene", ]
message("[generate_saf] ", nrow(gtf), " gene entries")

# ── Extract gene_id ──
extract_gene_id <- function(attr) {
    m <- regmatches(attr, regexec('gene_id "([^"]+)"', attr))[[1]]
    if (length(m) >= 2) return(m[2])
    m <- regmatches(attr, regexec('gene_id ([^;]+)', attr))[[1]]
    if (length(m) >= 2) return(trimws(m[2]))
    return(NA)
}

gene_ids <- sapply(gtf$attributes, extract_gene_id, USE.NAMES = FALSE)

# ── Build TSS SAF ──
tss_saf <- data.frame(
    GeneID = gene_ids, Chr = gtf$chr,
    Start  = integer(nrow(gtf)), End = integer(nrow(gtf)),
    Strand = gtf$strand, stringsAsFactors = FALSE
)

plus  <- gtf$strand == "+"
minus <- gtf$strand == "-"

tss_saf$Start[plus]  <- gtf$start[plus] - upstream
tss_saf$End[plus]    <- gtf$start[plus] + downstream
tss_saf$Start[minus] <- gtf$end[minus] - downstream
tss_saf$End[minus]   <- gtf$end[minus] + upstream

tss_saf$Start <- pmax(tss_saf$Start, 1)

# ── Build gene body SAF ──
gb_saf <- data.frame(
    GeneID = gene_ids, Chr = gtf$chr,
    Start  = integer(nrow(gtf)), End = integer(nrow(gtf)),
    Strand = gtf$strand, stringsAsFactors = FALSE
)

gb_saf$Start[plus]  <- gtf$start[plus] + offset
gb_saf$End[plus]    <- gtf$end[plus]
gb_saf$Start[minus] <- gtf$start[minus]
gb_saf$End[minus]   <- gtf$end[minus] - offset

gb_saf$Start <- pmax(gb_saf$Start, 1)

# ── Filter invalid ──
tss_ok <- !is.na(tss_saf$GeneID) & tss_saf$GeneID != "" & tss_saf$Start < tss_saf$End
gb_ok  <- !is.na(gb_saf$GeneID)  & gb_saf$GeneID  != "" & gb_saf$Start  < gb_saf$End

tss_saf <- tss_saf[tss_ok, ]
gb_saf  <- gb_saf[gb_ok, ]

# Deduplicate
tss_saf <- tss_saf[!duplicated(tss_saf$GeneID), ]
gb_saf  <- gb_saf[!duplicated(gb_saf$GeneID), ]

# Sort
chr_sort <- function(saf) {
    saf[order(suppressWarnings(as.integer(gsub("chr|CHR|Chr", "", saf$Chr))), saf$Chr, saf$Start), ]
}
tss_saf <- chr_sort(tss_saf)
gb_saf  <- chr_sort(gb_saf)

# ── Write outputs ──
write.table(tss_saf, file = argv$tss_out, sep = "\t", quote = FALSE, row.names = FALSE)
message("[generate_saf] TSS SAF: ", nrow(tss_saf), " regions → ", argv$tss_out)

write.table(gb_saf, file = argv$gb_out, sep = "\t", quote = FALSE, row.names = FALSE)
message("[generate_saf] Gene body SAF: ", nrow(gb_saf), " regions → ", argv$gb_out)

message("[generate_saf] Done.")
