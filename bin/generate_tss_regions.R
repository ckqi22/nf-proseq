#!/usr/bin/env Rscript
# =============================================================================
# generate_tss_regions.R — Generate TSS window SAF annotation from GTF
#
# Creates an SAF-format annotation file with a window around each gene's TSS.
# SAF format (used by featureCounts -F SAF): GeneID, Chr, Start, End, Strand
#
# For plus-strand genes:  window = (TSS - upstream) to (TSS + downstream)
# For minus-strand genes: window = (TSS - downstream) to (TSS + upstream)
# =============================================================================
suppressWarnings(suppressMessages({
    library(argparser)
}))

argv <- arg_parser("Generate TSS window SAF annotation from GTF for featureCounts")
argv <- add_argument(argv, "--gtf",        help = "Path to reference GTF file")
argv <- add_argument(argv, "--upstream",   help = "Upstream window size (bp)",   default = 50)
argv <- add_argument(argv, "--downstream", help = "Downstream window size (bp)", default = 300)
argv <- add_argument(argv, "--output",     help = "Output SAF file path")
argv <- parse_args(argv)

upstream   <- as.integer(argv$upstream)
downstream <- as.integer(argv$downstream)

if (!file.exists(argv$gtf)) stop("GTF file not found: ", argv$gtf)

message("[generate_tss_regions] Reading GTF: ", argv$gtf)
message("[generate_tss_regions] Window: -", upstream, " / +", downstream, " bp around TSS")

# Read GTF as a table (skip comment lines)
gtf <- read.delim(argv$gtf, header = FALSE, stringsAsFactors = FALSE, comment.char = "#")
names(gtf) <- c("chr", "source", "feature", "start", "end", "score", "strand", "frame", "attributes")

# Filter to transcripts
gtf <- gtf[gtf$feature %in% c("transcript", "mRNA"), ]
if (nrow(gtf) == 0) {
    # Try gene level if no transcripts
    gtf <- read.delim(argv$gtf, header = FALSE, stringsAsFactors = FALSE, comment.char = "#")
    names(gtf) <- c("chr", "source", "feature", "start", "end", "score", "strand", "frame", "attributes")
    gtf <- gtf[gtf$feature == "gene", ]
}
message("[generate_tss_regions] Found ", nrow(gtf), " transcript/gene entries")

# Extract gene_id from attributes column
extract_gene_id <- function(attr) {
    m <- regmatches(attr, regexec('gene_id "([^"]+)"', attr))[[1]]
    if (length(m) >= 2) return(m[2])
    m <- regmatches(attr, regexec('gene_id ([^;]+)', attr))[[1]]
    if (length(m) >= 2) return(trimws(m[2]))
    return(NA)
}

gene_ids <- sapply(gtf$attributes, extract_gene_id, USE.NAMES = FALSE)

# Determine TSS and create SAF window
saf <- data.frame(
    GeneID = gene_ids,
    Chr    = gtf$chr,
    Start  = integer(nrow(gtf)),
    End    = integer(nrow(gtf)),
    Strand = gtf$strand,
    stringsAsFactors = FALSE
)

# Plus strand: TSS = start
plus_idx <- gtf$strand == "+"
saf$Start[plus_idx] <- gtf$start[plus_idx] - upstream
saf$End[plus_idx]   <- gtf$start[plus_idx] + downstream

# Minus strand: TSS = end
minus_idx <- gtf$strand == "-"
saf$Start[minus_idx] <- gtf$end[minus_idx] - downstream
saf$End[minus_idx]   <- gtf$end[minus_idx] + upstream

# Clamp to positive coordinates
saf$Start <- pmax(saf$Start, 1)

# Remove entries with invalid gene IDs or coordinates
valid <- !is.na(saf$GeneID) & saf$GeneID != "" & saf$Start < saf$End
saf <- saf[valid, ]

# Deduplicate by GeneID (keep first occurrence = longest transcript usually)
saf <- saf[!duplicated(saf$GeneID), ]

# Sort by chromosome and start position
chr_order <- order(suppressWarnings(as.integer(gsub("chr|CHR|Chr", "", saf$Chr))),
                   saf$Chr, saf$Start)
saf <- saf[chr_order, ]

message("[generate_tss_regions] ", nrow(saf), " unique TSS regions generated")
message("[generate_tss_regions] Writing SAF: ", argv$output)
write.table(saf, file = argv$output, sep = "\t", quote = FALSE, row.names = FALSE)
message("[generate_tss_regions] Done.")
