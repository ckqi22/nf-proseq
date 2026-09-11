#!/usr/bin/env Rscript
# =============================================================================
# pausing_xlsx.R — 把 Pausing_Index.tsv + promoter.bed/genebody.bed + 注释表
#                   打包成 PROSeq_pausing.xlsx
#
# 每个样本一个 sheet（sheet 名 = 样本名），所有 sheet 列名与列序完全一致：
#   GeneID, chr, start, end, strand, Length, Annotation,
#   Pausing Index, Promoter Reads(-{upstream},{downstream}), GeneBody Reads({genebody_offset},TES)
#
# 列来源：
#   GeneID     = Pausing_Index.tsv 的 gene_id
#   chr/strand = promoter.bed（代表 transcript，第 1/6 列）
#   start      = TSS - upstream（promoter 5' 边界，取自 promoter.bed）
#   end        = TES（基因 3' 端，取自 genebody.bed）
#   Length     = end - start + 1（1-based 闭区间碱基数）
#   Annotation = 注释表第 2 列（gene_name）；无注释表时留空
#   Pausing Index            = {sample}_PI
#   Promoter Reads(-U,D)     = {sample}_promoter（U/D 随 --tss_upstream/--tss_downstream）
#   GeneBody Reads(O,TES)    = {sample}_genebody（O 随 --genebody_offset）
#
# 坐标约定（与 gtf2bed.R 一致）：promoter.bed/genebody.bed 为 BED6 0-based 半开，
#   且已按 strand 正确给出每基因的 promoter / genebody 窗口（基因组升序 Start<End）。
#   - promoter.bed = [TSS-upstream, TSS+downstream]（转录方向）
#   - genebody.bed = [TSS+genebody_offset, TES]（转录方向）
#   promoter 5' 边界（TSS-upstream）正链取 promoter.bed Start+1、负链取 End；
#   TES 正链取 genebody.bed End、负链取 Start+1。
#
# Usage:
#   Rscript pausing_xlsx.R --pi Pausing_Index.tsv \
#       --promoter_bed promoter.bed --genebody_bed genebody.bed \
#       [--annotation gene_annotation.txt] \
#       [--tss_upstream 50] [--tss_downstream 300] [--genebody_offset 301] \
#       --output PROSeq_pausing.xlsx
# =============================================================================
suppressPackageStartupMessages({
    library(argparse)
    library(openxlsx)
})

parser <- ArgumentParser(description = "Build PROSeq_pausing.xlsx (per-sample sheets)")
parser$add_argument("--pi",              required = TRUE,  help = "Pausing_Index.tsv")
parser$add_argument("--promoter_bed",    required = TRUE,  help = "promoter.bed (BED6, no header)")
parser$add_argument("--genebody_bed",    required = TRUE,  help = "genebody.bed (BED6, no header)")
parser$add_argument("--annotation",                        help = "Gene annotation table (first col gene_id, second gene_name). Optional.")
parser$add_argument("--tss_upstream",    type = "integer", default = 50,  help = "TSS upstream window (bp)")
parser$add_argument("--tss_downstream",  type = "integer", default = 300, help = "TSS downstream window (bp)")
parser$add_argument("--genebody_offset", type = "integer", default = 301, help = "Gene body start offset (bp)")
parser$add_argument("--output",          required = TRUE,  help = "Output xlsx path")
args <- parser$parse_args()

upstream        <- args$tss_upstream
downstream      <- args$tss_downstream
genebody_offset <- args$genebody_offset

# ── 读 Pausing_Index.tsv ──
pi <- read.delim(args$pi, header = TRUE, sep = "\t",
                 stringsAsFactors = FALSE, check.names = FALSE)

# ── 读 BED6（无表头，0-based 半开）──
read_bed <- function(path) {
    b <- read.delim(path, header = FALSE, sep = "\t",
                    stringsAsFactors = FALSE, comment.char = "#")
    colnames(b) <- c("chr", "start", "end", "GeneID", "score", "strand")
    b[!duplicated(b$GeneID), , drop = FALSE]   # gene 级唯一（GeneID 去重保首个）
}
prom_bed <- read_bed(args$promoter_bed)
gb_bed   <- read_bed(args$genebody_bed)

# ── 读注释表（可选）──
gene_name <- NULL
if (!is.null(args$annotation) && nzchar(args$annotation) && file.exists(args$annotation)) {
    ann <- read.delim(args$annotation, header = TRUE, sep = "\t",
                      stringsAsFactors = FALSE, check.names = FALSE, quote = "", comment.char = "")
    ann[is.na(ann)] <- ""
    # 优先 gene_name 列；否则取第 2 列
    name_col <- if ("gene_name" %in% colnames(ann)) "gene_name" else colnames(ann)[2]
    if (ncol(ann) >= 2) {
        gene_name <- setNames(as.character(ann[[name_col]]), as.character(ann[[1]]))
    }
}

# ── 样本名 = *_PI 列去掉后缀 ──
pi_cols <- grep("_PI$", colnames(pi), value = TRUE)
if (length(pi_cols) == 0)
    stop("[pausing_xlsx] no *_PI column found in ", args$pi)
samples <- sub("_PI$", "", pi_cols)

# ── 坐标 / 注释（按 gene_id 对齐到代表 transcript 的 promoter/genebody）──
pidx <- match(pi$gene_id, prom_bed$GeneID)
gidx <- match(pi$gene_id, gb_bed$GeneID)
if (anyNA(pidx) || anyNA(gidx))
    stop("[pausing_xlsx] some gene_id in Pausing_Index.tsv missing from promoter/genebody BED")

chr    <- prom_bed$chr[pidx]
strand <- prom_bed$strand[pidx]
plus   <- strand == "+"

# TSS-upstream（promoter 5' 边界，1-based）：
#   正链 = promoter.bed Start+1（左边界）；负链 = promoter.bed End（右边界）
prom_5p <- ifelse(plus, prom_bed$start[pidx] + 1, prom_bed$end[pidx])
# TES（基因 3' 端，1-based）：
#   正链 = genebody.bed End（右边界）；负链 = genebody.bed Start+1（左边界）
tes <- ifelse(plus, gb_bed$end[gidx], gb_bed$start[gidx] + 1)

# start/end 取基因组升序（保证 start <= end）；Length = end - start + 1（1-based 闭区间）
start <- pmin(prom_5p, tes)
end   <- pmax(prom_5p, tes)
span  <- end - start + 1

if (is.null(gene_name)) {
    annotation <- rep("", nrow(pi))
} else {
    annotation <- gene_name[pi$gene_id]
    annotation[is.na(annotation)] <- ""
}

# 动态列名（随窗口参数变化）
prom_col <- paste0("Promoter Reads(-", upstream, ",", downstream, ")")
gb_col   <- paste0("GeneBody Reads(", genebody_offset, ",TES)")

# ── 逐样本写 sheet ──
wb <- createWorkbook()
for (s in samples) {
    df <- data.frame(
        GeneID       = pi$gene_id,
        chr          = chr,
        start        = start,
        end          = end,
        strand       = strand,
        Length       = span,
        Annotation   = annotation,
        check.names  = FALSE,
        stringsAsFactors = FALSE
    )
    df[["Pausing Index"]] <- pi[[paste0(s, "_PI")]]
    df[[prom_col]]        <- pi[[paste0(s, "_promoter")]]
    df[[gb_col]]          <- pi[[paste0(s, "_genebody")]]

    # Excel sheet 名去非法字符 + 截断 31
    sname <- substr(gsub("[\\[\\]:*?/\\\\]", "_", s), 1, 31)
    addWorksheet(wb, sname)
    writeData(wb, sname, df, rowNames = FALSE)
    # 加粗表头 + 冻结首行
    addStyle(wb, sname, createStyle(textDecoration = "bold", halign = "center"),
             rows = 1, cols = 1:ncol(df))
    freezePane(wb, sname, firstActiveRow = 2)
    message("[pausing_xlsx] sheet '", sname, "': ", nrow(df), " genes")
}

saveWorkbook(wb, args$output, overwrite = TRUE)
message("[pausing_xlsx] wrote ", args$output, " (", length(samples), " sheets)")
