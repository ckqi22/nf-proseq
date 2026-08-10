#!/usr/bin/env Rscript
# =============================================================================
# gsea.R — Gene Set Enrichment Analysis using clusterProfiler
#
# Reads count matrix, computes log2FC per comparison, runs pre-ranked GSEA.
# =============================================================================
suppressWarnings(suppressMessages({
    library(argparser)
    library(clusterProfiler)
    library(org.Hs.eg.db)
    library(edgeR)
    library(openxlsx)
    library(ggplot2)
    library(enrichplot)
}))

argv <- arg_parser("GSEA analysis for PRO-seq data")
argv <- add_argument(argv, "--counts",      help = "Gene body count matrix (tab-separated)")
argv <- add_argument(argv, "--groups",      help = "YAML-like file mapping group names to sample lists")
argv <- add_argument(argv, "--comparisons", help = "CSV file with columns: group1, group2")
argv <- add_argument(argv, "--org_db",      help = "OrgDb package name", default = "org.Hs.eg.db")
argv <- add_argument(argv, "--organism",    help = "KEGG organism code", default = "hsa")
argv <- add_argument(argv, "--output_dir",  help = "Output directory for GSEA results")
argv <- parse_args(argv)

dir.create(argv$output_dir, recursive = TRUE, showWarnings = FALSE)

# --- Read count matrix ---
message("[gsea] Reading count matrix: ", argv$counts)
counts <- read.delim(argv$counts, header = TRUE, row.names = 1, stringsAsFactors = FALSE)
message("  ", nrow(counts), " genes x ", ncol(counts), " samples")

# --- Read groups ---
message("[gsea] Reading groups: ", argv$groups)
groups_raw <- read.delim(argv$groups, header = FALSE, stringsAsFactors = FALSE,
                         comment.char = "#", sep = ":", strip.white = TRUE)
group_names  <- trimws(groups_raw[, 1])
group_samples <- lapply(strsplit(trimws(groups_raw[, 2]), ","), trimws)
names(group_samples) <- group_names

sample_groups <- character()
for (gp in names(group_samples)) {
    for (samp in group_samples[[gp]]) {
        sample_groups[samp] <- gp
    }
}

# Match to count columns
count_samples_clean <- basename(gsub("\\.bam$", "", colnames(counts)))
group_vec <- factor(sample_groups[count_samples_clean])
valid <- !is.na(group_vec)
counts <- counts[, valid, drop = FALSE]
group_vec <- group_vec[valid]

# --- Read comparisons ---
message("[gsea] Reading comparisons: ", argv$comparisons)
comps <- read.delim(argv$comparisons, header = TRUE, stringsAsFactors = FALSE, sep = ",")

# --- Load OrgDb ---
org_db <- get(argv$org_db, envir = asNamespace(argv$org_db))

# --- Run GSEA for each comparison ---
for (i in seq_len(nrow(comps))) {
    g1 <- comps$group1[i]
    g2 <- comps$group2[i]
    comp_name <- paste0(g1, "_vs_", g2)

    message("[gsea] Processing: ", comp_name)

    if (!g1 %in% levels(group_vec) || !g2 %in% levels(group_vec)) {
        message("  Skipping: group not found")
        next
    }

    # Subset to the two groups
    idx <- group_vec %in% c(g1, g2)
    sub_counts <- counts[, idx, drop = FALSE]
    sub_group  <- factor(group_vec[idx], levels = c(g2, g1))  # g2 = reference

    # edgeR DGEList
    dge <- DGEList(counts = sub_counts, group = sub_group)
    keep <- filterByExpr(dge, group = sub_group)
    dge <- dge[keep, , keep.lib.sizes = FALSE]
    dge <- calcNormFactors(dge, method = "TMM")

    # Estimate dispersion and run exact test
    design <- model.matrix(~ sub_group)
    dge <- estimateDisp(dge, design)
    et <- exactTest(dge, pair = c(g2, g1))

    res <- topTags(et, n = Inf)$table
    res$gene_id <- rownames(res)

    # Convert to ENTREZ and rank by logFC
    gene_map <- tryCatch(
        bitr(res$gene_id, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org_db),
        error = function(e) NULL
    )

    if (is.null(gene_map) || nrow(gene_map) == 0) {
        message("  No valid ENTREZ IDs; skipping")
        next
    }

    res_mapped <- merge(res, gene_map, by.x = "gene_id", by.y = "SYMBOL")

    # Rank by logFC (descending)
    gene_list <- res_mapped$logFC
    names(gene_list) <- res_mapped$ENTREZID
    gene_list <- sort(gene_list, decreasing = TRUE)

    # GO GSEA
    for (ont in c("BP", "MF", "CC")) {
        gsea_go <- tryCatch(
            gseGO(
                geneList     = gene_list,
                ont          = ont,
                OrgDb        = org_db,
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                verbose       = FALSE
            ),
            error = function(e) { message("  WARNING: gseGO (", ont, "): ", e$message); NULL }
        )

        if (!is.null(gsea_go) && nrow(gsea_go) > 0) {
            prefix <- paste0(comp_name, "_GO_", ont)
            write.xlsx(as.data.frame(gsea_go),
                       file = file.path(argv$output_dir, paste0(prefix, ".xlsx")),
                       rowNames = FALSE)

            # Enrichment plot
            pdf(file.path(argv$output_dir, paste0(prefix, "_gseaplot.pdf")), width = 10, height = 8)
            top_pathways <- min(5, nrow(gsea_go))
            for (j in seq_len(top_pathways)) {
                print(gseaplot2(gsea_go, geneSetID = j, title = gsea_go$Description[j]))
            }
            dev.off()

            message("  GO_", ont, ": ", nrow(gsea_go), " enriched gene sets")
        }
    }

    # KEGG GSEA
    gsea_kegg <- tryCatch(
        gseKEGG(
            geneList     = gene_list,
            organism     = argv$organism,
            pAdjustMethod = "BH",
            pvalueCutoff  = 0.05,
            verbose       = FALSE
        ),
        error = function(e) { message("  WARNING: gseKEGG: ", e$message); NULL }
    )

    if (!is.null(gsea_kegg) && nrow(gsea_kegg) > 0) {
        prefix <- paste0(comp_name, "_KEGG")
        write.xlsx(as.data.frame(gsea_kegg),
                   file = file.path(argv$output_dir, paste0(prefix, ".xlsx")),
                   rowNames = FALSE)

        pdf(file.path(argv$output_dir, paste0(prefix, "_gseaplot.pdf")), width = 10, height = 8)
        top_pathways <- min(5, nrow(gsea_kegg))
        for (j in seq_len(top_pathways)) {
            print(gseaplot2(gsea_kegg, geneSetID = j, title = gsea_kegg$Description[j]))
        }
        dev.off()

        message("  KEGG: ", nrow(gsea_kegg), " enriched pathways")
    }
}

message("[gsea] Done.")
