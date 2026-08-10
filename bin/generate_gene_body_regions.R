#!/usr/bin/env Rscript
# =============================================================================
# generate_gene_body_regions.R — Generate gene body SAF from GTF
#
# Gene body = TSS+offset to TES (excludes TSS-proximal pause region).
# SAF format (featureCounts -F SAF): GeneID, Chr, Start, End, Strand
#
# Plus-strand genes:  TSS = start,  TES = end   → body = (start+offset) to end
# Minus-strand genes: TSS = end,    TES = start → body = start to (end-offset)
# =============================================================================
suppressWarnings(suppressMessages({
    library(argparser)
}))

argv <- arg_parser("Generate gene body SAF annotation from GTF for featureCounts")
argv <- add_argument(argv, "--gtf",    help = "Path to reference GTF file")
argv <- add_argument(argv, "--offset", help = "Distance from TSS to gene body start (bp)", default = 301)
argv <- add_argument(argv, "--output", help = "Output SAF file path")
argv <- parse_args(argv)

offset <- as.integer(argv$offset)

if (!file.exists(argv$gtf)) stop("GTF file not found: ", argv$gtf)

message("[generate_gene_body] Reading GTF: ", argv$gtf)
message("[generate_gene_body] Gene body offset: TSS+", offset, "bp")

# Read GTF
gtf <- read.delim(argv$gtf, header = FALSE, stringsAsFactors = FALSE, comment.char = "#")
names(gtf) <- c("chr", "source", "feature", "start", "end", "score", "strand", "frame", "attributes")

# Filter to transcripts; fall back to genes
gtf <- gtf[gtf$feature %in% c("transcript", "mRNA"), ]
if (nrow(gtf) == 0) {
    gtf <- read.delim(argv$gtf, header = FALSE, stringsAsFactors = FALSE, comment.char = "#")
    names(gtf) <- c("chr", "source", "feature", "start", "end", "score", "strand", "frame", "attributes")
    gtf <- gtf[gtf$feature == "gene", ]
}
message("[generate_gene_body] Found ", nrow(gtf), " transcript/gene entries")

# Extract gene_id
extract_gene_id <- function(attr) {
    m <- regmatches(attr, regexec('gene_id "([^"]+)"', attr))[[1]]
    if (length(m) >= 2) return(m[2])
    m <- regmatches(attr, regexec('gene_id ([^;]+)', attr))[[1]]
    if (length(m) >= 2) return(trimws(m[2]))
    return(NA)
}

gene_ids <- sapply(gtf$attributes, extract_gene_id, USE.NAMES = FALSE)

# Build SAF: TSS+offset to TES
saf <- data.frame(
    GeneID = gene_ids,
    Chr    = gtf$chr,
    Start  = integer(nrow(gtf)),
    End    = integer(nrow(gtf)),
    Strand = gtf$strand,
    stringsAsFactors = FALSE
)

# Plus strand: TSS = start, TES = end
plus_idx <- gtf$strand == "+"
saf$Start[plus_idx] <- gtf$start[plus_idx] + offset
saf$End[plus_idx]   <- gtf$end[plus_idx]

# Minus strand: TSS = end, TES = start
minus_idx <- gtf$strand == "-"
saf$Start[minus_idx] <- gtf$start[minus_idx]
saf$End[minus_idx]   <- gtf$end[minus_idx] - offset

# Clamp: start must be >= 1
saf$Start <- pmax(saf$Start, 1)

# Remove invalid: gene_id missing, or start >= end (gene too short)
valid <- !is.na(saf$GeneID) & saf$GeneID != "" & saf$Start < saf$End
n_removed <- sum(!valid)
if (n_removed > 0) {
    message("[generate_gene_body] Removed ", n_removed, " entries (invalid or too short for offset=", offset, ")")
}
saf <- saf[valid, ]

# Deduplicate by GeneID (keep first = longest transcript)
saf <- saf[!duplicated(saf$GeneID), ]

# Sort
chr_order <- order(suppressWarnings(as.integer(gsub("chr|CHR|Chr", "", saf$Chr))),
                   saf$Chr, saf$Start)
saf <- saf[chr_order, ]

message("[generate_gene_body] ", nrow(saf), " gene body regions generated")
message("[generate_gene_body] Writing SAF: ", argv$output)
write.table(saf, file = argv$output, sep = "\t", quote = FALSE, row.names = FALSE)
message("[generate_gene_body] Done.")
