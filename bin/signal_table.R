#!/usr/bin/env Rscript
# =============================================================================
# signal_table.R — 活跃蛋白编码基因 pause 窗口的单碱基 Pol II 信号表
#
# 筛选级联（阈值自适应测序深度/标准化方式，无绝对信号常数）：
#   基因级：蛋白编码 + 代表转录本（longest_tx.R 已完成）→ 活跃（任一组）：
#           promoter 密度 > 0 且 genebody 密度 >= 候选基因密度中位数 × active_frac
#   区域级：pause 窗口 = promoter.bed（gtf2bed.R 产物，params.yml 单一口径）
#   位点级：组均值 Signal >= max(本基因峰值 × peak_frac, 噪声线)，
#           噪声线 = 全体窗口位点 Signal 的 noise_quantile 分位；
#           且组内 >= min_reps 个重复达同阈值（重复不足时要求全部）
#   位点归属：divergent promoter 按 (位点, 基因) 对输出；组间 union 对齐，缺失填 0
# =============================================================================
suppressWarnings(suppressMessages({
    library(argparser)
    library(rtracklayer)
    library(GenomicRanges)
}))

# ── 参数（全部 CLI 外放）──
argv <- arg_parser("Per-base Pol II signal table of active genes' pause windows")
argv <- add_argument(argv, "--manifest",       help = "TSV no header, 6 cols: group \\t sample \\t plus_raw \\t minus_raw \\t plus_cpm \\t minus_cpm (bigWig paths)")
argv <- add_argument(argv, "--promoter_bed",   help = "BED6 pause windows (gtf2bed.R promoter.bed, 0-based half-open)")
argv <- add_argument(argv, "--promoter_count", help = "pol2_promoter.matrix.txt (gene_id, length, samples...)")
argv <- add_argument(argv, "--genebody_count", help = "pol2_genebody.matrix.txt (gene_id, length, samples...)")
argv <- add_argument(argv, "--rep_gtf",        help = "longest_tx.gtf (representative transcripts) -> transcriptid column")
argv <- add_argument(argv, "--annotation",     help = "Annotation table, first col gene_id; remaining cols appended to the table")
argv <- add_argument(argv, "--active_frac",    help = "Active-gene genebody density floor = median density of candidate genes x this fraction (default: 0.005 ~ 0.2 CPM/kb on fig4, ~ NRSA 4 reads/kb at 20M lib)", default = 0.005)
argv <- add_argument(argv, "--min_genebody_length", help = "Min genebody region length, bp (default: 800, same as pausing_index.R)", default = 800)
argv <- add_argument(argv, "--peak_frac",      help = "Site threshold: group-mean Signal >= gene pause-peak x this fraction (default: 0.1)", default = 0.1)
argv <- add_argument(argv, "--noise_quantile", help = "Site noise floor = this quantile of all pause-window site Signals (default: 0.9; 0 = off)", default = 0.9)
argv <- add_argument(argv, "--min_reps",       help = "Replicate support: >= this many replicates with Signal >= site threshold (default: 2; all required if group has fewer)", default = 2)
argv <- add_argument(argv, "--orientation",    help = "Library orientation: 'reverse' (R1 antisense, active site = read 5' end) or 'forward' (R1 sense, active site = read 3' end)", default = "reverse")
argv <- add_argument(argv, "--output",         help = "Output directory (writes pol2_signal_table.tsv + pol2_signal_table.note.txt) (default: ./)", default = "./")
argv <- parse_args(argv)

# ── 输入读取 ──

read_manifest <- function(path) {
    man <- read.delim(path, header = FALSE, sep = "\t", stringsAsFactors = FALSE,
                      col.names = c("group", "sample", "plus_raw", "minus_raw",
                                    "plus_cpm", "minus_cpm"))
    if (nrow(man) == 0) stop("[signal_table] empty manifest")
    man
}

read_promoter_bed <- function(path, valid_chrom) {
    g <- read.delim(path, header = FALSE, sep = "\t", stringsAsFactors = FALSE,
                    col.names = c("chrom", "start", "end", "gene_id", "score", "strand"))
    if (any(duplicated(g$gene_id))) {
        warning("[signal_table] duplicated gene_id in promoter.bed, keeping first", call. = FALSE)
        g <- g[!duplicated(g$gene_id), , drop = FALSE]
    }
    g[g$chrom %in% valid_chrom, , drop = FALSE]
}

read_count_matrix <- function(path, label) {
    m <- read.delim(path, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
                    check.names = FALSE, quote = "")
    if (!all(c("gene_id", "length") %in% names(m)))
        stop("[signal_table] ", label, " must have gene_id and length columns")
    list(gene_id = m$gene_id, length = as.numeric(m$length),
         counts = as.matrix(m[, setdiff(names(m), c("gene_id", "length")), drop = FALSE]))
}

# 两矩阵按 gene_id 对齐（length 列本就不同：promoter 窗口恒定 vs genebody 可变）
align_count_matrices <- function(prom, gbd) {
    if (identical(prom$gene_id, gbd$gene_id)) return(list(prom = prom, gbd = gbd))
    common <- intersect(prom$gene_id, gbd$gene_id)
    warning("[signal_table] promoter/genebody matrices differ: keeping ",
            length(common), " common genes", call. = FALSE)
    i1 <- match(common, prom$gene_id); i2 <- match(common, gbd$gene_id)
    list(prom = list(gene_id = common, length = prom$length[i1], counts = prom$counts[i1, , drop = FALSE]),
         gbd  = list(gene_id = common, length = gbd$length[i2],  counts = gbd$counts[i2, , drop = FALSE]))
}

# ── 基因级：密度 → 活跃判定 ──

# 密度 = count / 文库代理 × 1e6 / kb（文库代理 = genebody 矩阵列和）
count_to_density <- function(counts, libsize, region_len_bp) {
    sweep(counts, 2, libsize, "/") * 1e6 / (region_len_bp / 1000)
}

group_mean <- function(dens, man, groups) {
    out <- list()
    for (grp in groups)
        out[[grp]] <- rowMeans(dens[, man$sample[man$group == grp], drop = FALSE])
    out
}

call_active_genes <- function(prom_g, gb_g, groups, active_frac) {
    active <- rep(FALSE, length(prom_g[[1]]))
    for (grp in groups) {
        thr_g <- median(gb_g[[grp]]) * active_frac
        active <- active | (prom_g[[grp]] > 0 & gb_g[[grp]] >= thr_g)
        message("[signal_table] group ", grp, ": genebody density median = ",
                signif(median(gb_g[[grp]]), 3), ", active threshold = ",
                signif(thr_g, 3), ", passing = ", sum(prom_g[[grp]] > 0 & gb_g[[grp]] >= thr_g))
    }
    message("[signal_table] active genes (any group): ", sum(active), " / ", length(active))
    active
}

# ── 位点级：单碱基信号提取与筛选 ──

# 单样本 pause 窗口内单碱基信号（raw + CPM 双轨，展开为碱基级）
# win_gr: 候选基因的 pause 窗口（仅含与信号文件同向的基因）
pause_window_bases <- function(bw_raw, bw_cpm, win_gr, win_gene_idx, label) {
    sig_raw <- rtracklayer::import(BigWigFile(bw_raw), which = win_gr)
    if (!length(sig_raw)) {
        warning("[signal_table] no signal found for ", label, call. = FALSE)
        return(NULL)
    }
    extract <- function(sig_gr) {
        hits <- findOverlaps(win_gr, sig_gr)
        if (!length(hits)) return(NULL)
        ir <- pintersect(win_gr[queryHits(hits)], sig_gr[subjectHits(hits)])
        lens <- width(ir)
        data.frame(chrom    = rep(as.character(seqnames(ir)), lens),
                   pos      = rep(start(ir), lens) + sequence(lens) - 1L,
                   gene_idx = rep(win_gene_idx[queryHits(hits)], lens),
                   value    = abs(rep(sig_gr$score[subjectHits(hits)], lens)),  # minus 轨存负值
                   stringsAsFactors = FALSE)
    }
    raw_df <- extract(sig_raw)
    if (is.null(raw_df)) return(NULL)
    cpm_df <- extract(rtracklayer::import(BigWigFile(bw_cpm), which = win_gr))
    key <- paste(raw_df$chrom, raw_df$pos, raw_df$gene_idx, sep = "|")
    if (is.null(cpm_df)) raw_df$value_cpm <- 0 else {
        cpm_key <- paste(cpm_df$chrom, cpm_df$pos, cpm_df$gene_idx, sep = "|")
        raw_df$value_cpm <- cpm_df$value[match(key, cpm_key)]
        raw_df$value_cpm[is.na(raw_df$value_cpm)] <- 0
    }
    raw_df$value_raw <- raw_df$value
    raw_df[, c("chrom", "pos", "gene_idx", "value_raw", "value_cpm")]
}

# 位点 universe：活跃基因 pause 窗口内全部样本单碱基并集（raw + CPM 双轨，缺失填 0）
build_site_universe <- function(man, g, active) {
    win_idx <- which(active)
    win_gr  <- GRanges(g$chrom[win_idx],
                       IRanges(pmax(g$start[win_idx], 0L) + 1L, g$end[win_idx]),
                       strand = g$strand[win_idx])
    sample_keys <- paste0(man$group, ".", man$sample)
    site_list <- list()
    for (i in seq_len(nrow(man))) {
        key <- sample_keys[i]
        for (s in c("plus", "minus")) {
            want <- ifelse(s == "plus", "+", "-")
            idx_s <- win_idx[g$strand[win_idx] == want]
            if (!length(idx_s)) next
            df <- pause_window_bases(man[[paste0(s, "_raw")]][i], man[[paste0(s, "_cpm")]][i],
                                     win_gr[g$strand[win_idx] == want], idx_s, paste0(key, " ", s))
            if (!is.null(df)) site_list[[paste0(key, ".", s)]] <- df
        }
    }
    all_sites <- do.call(rbind, site_list)
    sites <- unique(all_sites[, c("chrom", "pos", "gene_idx")])
    sites$strand <- g$strand[sites$gene_idx]
    message("[signal_table] site universe (per-base union): ", nrow(sites))

    key_of_site <- paste(sites$chrom, sites$pos, sites$gene_idx, sep = "|")
    raw_mat <- matrix(0, nrow = nrow(sites), ncol = nrow(man), dimnames = list(NULL, sample_keys))
    cpm_mat <- raw_mat
    for (i in seq_len(nrow(man))) {
        key <- sample_keys[i]
        for (s in c("plus", "minus")) {
            df <- site_list[[paste0(key, ".", s)]]
            if (is.null(df)) next
            m <- match(paste(df$chrom, df$pos, df$gene_idx, sep = "|"), key_of_site)
            has <- !is.na(m)
            raw_mat[m[has], i] <- df$value_raw[has]
            cpm_mat[m[has], i] <- df$value_cpm[has]
        }
    }
    list(sites = sites, raw_mat = raw_mat, cpm_mat = cpm_mat)
}

# 位点阈值 = max(基因峰值 × peak_frac, 噪声线)；返回 list(thr, noise_floor)
site_thresholds <- function(mean_sig, sites, peak_frac, noise_q) {
    site_max  <- do.call(pmax, c(as.data.frame(as.matrix(mean_sig)), list(na.rm = TRUE)))
    gene_peak <- tapply(site_max, sites$gene_idx, max)          # 命名向量（名字 = gene_idx）
    noise_floor <- if (noise_q > 0) as.numeric(quantile(site_max, noise_q)) else 0
    message("[signal_table] noise floor (P", round(noise_q * 100),
            " of site Signals) = ", signif(noise_floor, 3))
    list(thr = pmax(as.numeric(gene_peak[as.character(sites$gene_idx)]) * peak_frac, noise_floor),
         noise_floor = noise_floor)
}

# 位点保留：任一组（组均值 >= 阈值 且 >= min_reps 重复 >= 阈值）
filter_sites <- function(sites, cpm_mat, man, groups, mean_sig, thr_site, min_reps) {
    keep <- rep(FALSE, nrow(sites))
    for (grp in groups) {
        cols <- which(man$group == grp)
        mean_pass <- mean_sig[[paste0(grp, ".Signal")]] >= thr_site
        n_pass <- rowSums(cpm_mat[, cols, drop = FALSE] >= thr_site)
        need <- min(min_reps, length(cols))   # 组内重复数不足时要求全部达标
        keep <- keep | (mean_pass & n_pass >= need)
    }
    message("[signal_table] sites passing threshold + replicate support: ",
            sum(keep), " / ", nrow(sites))
    keep
}

# ── 注释与输出 ──

# rep_gtf：gene_id → transcript_id（代表转录本）
load_transcript_map <- function(rep_gtf) {
    lines <- readLines(rep_gtf, warn = FALSE)
    lines <- lines[!grepl("^#", lines) & nzchar(lines)]
    attr_of <- function(line, key) {
        m <- regmatches(line, regexpr(paste0(key, ' "[^"]*"'), line))
        if (!length(m)) NA_character_ else sub(paste0(key, ' "'), "", substr(m, 1, nchar(m) - 1))
    }
    parts <- strsplit(lines, "\t", fixed = TRUE)
    tx_map <- unique(data.frame(
        gene_id = vapply(parts, function(p) attr_of(p[9], "gene_id"), character(1)),
        transcriptid = vapply(parts, function(p) attr_of(p[9], "transcript_id"), character(1)),
        stringsAsFactors = FALSE))
    tx_map <- tx_map[!is.na(tx_map$gene_id) & !is.na(tx_map$transcriptid), , drop = FALSE]
    tx_map[!duplicated(tx_map$gene_id), , drop = FALSE]
}

# 注释表连接（首列 gene_id，其余列透传表尾；NA 填空串）
append_annotation <- function(out, annotation_path) {
    ann <- read.delim(annotation_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
                      check.names = FALSE, quote = "")
    names(ann)[1] <- "gene_id"
    ann[is.na(ann)] <- ""
    ann <- ann[!duplicated(ann$gene_id), , drop = FALSE]
    if (!nrow(ann)) warning("[signal_table] empty annotation table", call. = FALSE)
    ai <- match(out$gene_id, ann$gene_id)
    for (cn in setdiff(names(ann), "gene_id")) {
        v <- ann[[cn]][ai]
        out[[cn]] <- ifelse(is.na(v), "", v)
    }
    out
}

# ── Main（薄编排层：读入 → 基因级 → 位点级 → 组装 → 输出）──
main <- function(argv) {
    man <- read_manifest(argv$manifest)
    groups <- unique(man$group)
    message("[signal_table] groups: ", paste(groups, collapse = ", "),
            " (", nrow(man), " samples)")

    # 基因级：promoter.bed ∩ 计数矩阵 → 密度 → 活跃
    si <- seqinfo(BigWigFile(man$plus_raw[1]))
    g <- read_promoter_bed(argv$promoter_bed, seqnames(si))
    mats <- align_count_matrices(read_count_matrix(argv$promoter_count, "promoter_count"),
                                 read_count_matrix(argv$genebody_count, "genebody_count"))
    if (!all(man$sample %in% colnames(mats$prom$counts)))
        stop("[signal_table] manifest samples (",
             paste(setdiff(man$sample, colnames(mats$prom$counts)), collapse = ", "),
             ") missing from count matrix columns (",
             paste(colnames(mats$prom$counts), collapse = ", "), ")")
    g <- g[g$gene_id %in% mats$prom$gene_id, , drop = FALSE]
    mi <- match(g$gene_id, mats$prom$gene_id)
    g$prom_len <- mats$prom$length[mi]; g$gb_len <- mats$gbd$length[mi]
    g <- g[g$gb_len >= as.integer(argv$min_genebody_length), , drop = FALSE]
    message("[signal_table] candidate genes (promoter.bed ∩ matrices, body >= ",
            argv$min_genebody_length, " bp): ", nrow(g))
    if (!nrow(g)) stop("[signal_table] no candidate genes")
    mi <- match(g$gene_id, mats$prom$gene_id)
    libsize   <- colSums(mats$gbd$counts)
    prom_dens <- count_to_density(mats$prom$counts[mi, man$sample, drop = FALSE], libsize, g$prom_len)
    gb_dens   <- count_to_density(mats$gbd$counts[mi, man$sample, drop = FALSE], libsize, g$gb_len)
    active <- call_active_genes(group_mean(prom_dens, man, groups),
                                group_mean(gb_dens, man, groups),
                                groups, as.numeric(argv$active_frac))
    if (!sum(active)) stop("[signal_table] no active genes survive the gene-level cascade")

    # 位点级：单碱基 universe → 阈值（峰形 + 噪声线）→ 重复支持度
    uni <- build_site_universe(man, g, active)
    mean_sig <- lapply(groups, function(grp)
        round(rowMeans(uni$cpm_mat[, which(man$group == grp), drop = FALSE]), 2))
    names(mean_sig) <- paste0(groups, ".Signal")
    mean_sig <- as.data.frame(mean_sig, stringsAsFactors = FALSE)
    thr <- site_thresholds(mean_sig, uni$sites,
                           as.numeric(argv$peak_frac), as.numeric(argv$noise_quantile))
    keep <- filter_sites(uni$sites, uni$cpm_mat, man, groups, mean_sig, thr$thr,
                         as.integer(argv$min_reps))
    if (!sum(keep)) stop("[signal_table] no sites survive the site-level cascade")

    # 组装输出表（用户定义列序）
    sites <- uni$sites[keep, , drop = FALSE]
    gsel  <- g[sites$gene_idx, , drop = FALSE]
    out <- data.frame(chrom = sites$chrom, start = sites$pos - 1L, end = sites$pos,   # 0-based 半开
                      stringsAsFactors = FALSE)
    # 各组 Count 集中在前、Signal 集中在后（跨组同类型列相邻）
    for (grp in groups) {
        cols <- which(man$group == grp)
        out[[paste0(grp, ".Count")]] <- round(rowSums(uni$raw_mat[keep, cols, drop = FALSE]))
    }
    for (grp in groups) {
        out[[paste0(grp, ".Signal")]] <- mean_sig[[paste0(grp, ".Signal")]][keep]
    }
    tx_map <- load_transcript_map(argv$rep_gtf)
    out$transcriptid <- tx_map$transcriptid[match(gsel$gene_id, tx_map$gene_id)]
    out$transcriptid[is.na(out$transcriptid)] <- ""
    out$gene_id <- gsel$gene_id
    out$Strand  <- gsel$strand
    out <- append_annotation(out, argv$annotation)
    chrom_levels <- seqnames(si)[seqnames(si) %in% unique(out$chrom)]
    out <- out[order(factor(out$chrom, levels = chrom_levels), out$start, out$Strand), ]

    note <- paste0(
        "Per-base Pol II signal table of active protein-coding genes' pause windows. ",
        "Gene standard: protein-coding + representative transcript (longest_tx.R) -> active ",
        "(promoter density > 0 AND genebody density >= median density of candidate genes x ",
        argv$active_frac, " in any group; fig4-calibrated ~0.2 CPM/kb, ~ NRSA 4 reads/kb at 20M lib). ",
        "Region: pause window = promoter.bed (params.yml single source). ",
        "Site standard: group-mean Signal >= max(gene pause-peak x ", argv$peak_frac,
        ", noise floor = P", round(as.numeric(argv$noise_quantile) * 100),
        " of all pause-window site Signals",
        if (as.numeric(argv$noise_quantile) > 0) paste0(" = ", signif(thr$noise_floor, 3)) else "",
        ") AND >= ", argv$min_reps,
        " replicates with Signal >= the same threshold in any group (all replicates required ",
        "if group has fewer). Count = sum of raw ",
        ifelse(identical(argv$orientation, "forward"), "3'-end", "5'-end"),
        " counts within group (integer evidence; active-site end per strandedness=", argv$orientation, "); ",
        "Signal = group-mean normalized signal (per-sample normalized then averaged; missing = 0). ",
        "Single-base, 0-based half-open; signal strand = gene strand. ",
        "Divergent-promoter sites appear once per gene (site, gene) pair. ",
        "Density normalization: count / genebody-matrix column sum x 1e6 / kb (library-size proxy).")

    dir.create(argv$output, showWarnings = FALSE, recursive = TRUE)
    out_path <- file.path(argv$output, "pol2_signal_table.tsv")
    write.table(out, file = out_path, sep = "\t", quote = FALSE, row.names = FALSE)
    writeLines(note, file.path(argv$output, "pol2_signal_table.note.txt"))
    message("[signal_table] ", nrow(out), " signal rows, ", length(unique(out$gene_id)),
            " genes -> ", out_path)
}

main(argv)
