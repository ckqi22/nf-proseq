#!/usr/bin/env bash
# ============================================================
# TT-Seq FastQC 原图提取：从 FastQC zip 提取两张交付图
#   per_base_quality.png（碱基质量分布）+ adapter_content.png（接头）
# 分 raw/trimmed 两个目录，供交付报告
#   03.Data_QC/FastQC/{Raw,Trimmed} 使用。
# FastQC zip 内结构固定为 {输入基础名}_fastqc/Images/*.png，
#   故用 unzip -p "*/Images/{img}.png" 提取（不展开整个 zip）。
# 用法: bash extract_fastqc_images.sh <raw_dir> <trimmed_dir> <outdir>
#   raw_dir     : fastqc/（fastp 前，原始 reads 的 FastQC zip）
#   trimmed_dir : fastqc_trimmed/（fastp 后，trimmed reads 的 FastQC zip）
#   outdir      : 输出根目录（其下生成 raw/ 与 trimmed/ 两个子目录）
# 输出命名: {sample}_R{r}_{per_base_quality|adapter_content}.png
# ============================================================
set -euo pipefail

if [ "$#" -ne 3 ]; then
    echo "用法: bash extract_fastqc_images.sh <raw_dir> <trimmed_dir> <outdir>" >&2
    exit 1
fi

RAW_DIR="$1"
TRIMMED_DIR="$2"
OUTDIR="$3"

extract_stage () {
    local zip_dir="$1" stage="$2"
    local dest="$OUTDIR/$stage"
    mkdir -p "$dest"
    local zip
    for zip in "$zip_dir"/*_fastqc.zip; do
        [ -e "$zip" ] || continue
        # 文件名 {sample}_R1_fastqc.zip / {sample}_R2_fastqc.zip → {sample}_R1 / {sample}_R2
        local base
        base=$(basename "$zip" .zip)          # {sample}_R1_fastqc
        local sample_r=${base%_fastqc}         # {sample}_R1
        local img
        for img in per_base_quality adapter_content; do
            unzip -p "$zip" "*/Images/${img}.png" > "$dest/${sample_r}_${img}.png" 2>/dev/null || true
        done
    done
}

extract_stage "$RAW_DIR"     raw
extract_stage "$TRIMMED_DIR" trimmed

echo "FastQC 原图提取完成:"
ls -1 "$OUTDIR"/raw "$OUTDIR"/trimmed 2>/dev/null || true
