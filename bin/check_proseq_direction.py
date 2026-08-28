#!/usr/bin/env python3
"""
确认 PRO-seq 建库的链方向与信号端（read 的 5' / 3' 端哪个是聚合酶位点）。

用途：新数据接入流程前，先用此脚本确认两个关键约定，再据此设置
  params.strandedness（read 相对转录方向）与 params.signal_end（聚合酶位点在
  read 的哪一端）。标准 PRO-seq small-RNA 库的预期结果是：
    read1 = antisense -> strandedness = reverse
    read1 5' 端 = 聚合酶位点 -> signal_end = 5

依赖：samtools（需先 source 到含 samtools 的 conda 环境）。

用法：
    python3 check_proseq_direction.py <BAM> <GTF> [--genes Actb,Gapdh,...]

输出：
    1) 每个基因 read1 比对到 + / - 链的计数 -> 判定 sense / antisense；
    2) 6 个基因合并的 TSS 相对位置 profile（下游为正）的 5' / 3' 端计数 ->
       依据暂停峰（TSS 下游 +20~+100bp）判定 signal_end。
"""

import argparse
import re
import subprocess
import sys
from collections import defaultdict

# 高表达看家基因（mouse mm39 已验证；其它物种用 --genes 覆盖）。
DEFAULT_GENES = ["Actb", "Gapdh", "Eef2", "Pgk1", "Rpl13a", "Ubc"]


def parse_gene(gtf, names):
    """从 GTF 提取每个基因的一条 gene 记录: name -> (chr, start, end, strand)."""
    wanted = set(names)
    found = {}
    with open(gtf) as f:
        for line in f:
            if line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 9 or c[2] != "gene":
                continue
            m = re.search(r'gene_name "([^"]+)"', c[8])
            if not m:
                continue
            g = m.group(1)
            if g in wanted and g not in found:
                found[g] = (c[0], int(c[3]), int(c[4]), c[6])
    return found


def cigar_len(c):
    n = 0
    for m in re.finditer(r"(\d+)([MIDNSHP=X])", c):
        if m.group(2) in "MDN=X":
            n += int(m.group(1))
    return n


def five_three(flag, pos, length):
    """按 FLAG 计算 read 的 5' 端与 3' 端基因组坐标（1-based）。"""
    if flag & 16:  # reverse strand
        return pos + length - 1, pos
    return pos, pos + length - 1


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("bam", help="比对后的 BAM（read1 单端或含 read1 的比对）")
    ap.add_argument("gtf", help="参考 GTF")
    ap.add_argument("--genes", default=",".join(DEFAULT_GENES),
                    help="用于判定的基因名（逗号分隔），须覆盖 + / - 两条链")
    args = ap.parse_args()

    genes = parse_gene(args.gtf, args.genes.split(","))
    missing = set(args.genes.split(",")) - set(genes)
    if missing:
        sys.exit(f"GTF 中未找到基因: {sorted(missing)}")

    # ---------- (1) 链方向 ----------
    print("========== (1) read1 比对链（sense / antisense） ==========")
    print(f"{'gene':8s} {'str':3s} | {'read1 +链':>9s} {'read1 -链':>9s} | 结论")
    n_anti = 0
    for name in args.genes.split(","):
        chr_, s, e, strand = genes[name]
        out = subprocess.check_output(
            ["samtools", "view", args.bam, f"{chr_}:{s}-{e}"], text=True)
        p = m = 0
        for line in out.splitlines():
            flag = int(line.split("\t")[1])
            if flag & 16:
                m += 1
            else:
                p += 1
        verdict = "sense" if ((strand == "+" and p > m) or
                              (strand == "-" and m > p)) else "antisense"
        if verdict == "antisense":
            n_anti += 1
        print(f"{name:8s} {strand:3s} | {p:9d} {m:9d} | {verdict}")

    n = len(args.genes.split(","))
    if n_anti == n:
        print(f"\n结论: read1 = antisense（{n_anti}/{n} 基因）→ strandedness = reverse")
    elif n_anti == 0:
        print(f"\n结论: read1 = sense（0/{n} antisense）→ strandedness = forward")
    else:
        print(f"\n警告: 结果不一致（{n_anti}/{n} antisense），请人工检查")

    # ---------- (2) 信号端：TSS profile ----------
    prof5 = defaultdict(int)
    prof3 = defaultdict(int)
    for name in args.genes.split(","):
        chr_, s, e, strand = genes[name]
        tss = s if strand == "+" else e
        out = subprocess.check_output(
            ["samtools", "view", args.bam, f"{chr_}:{tss - 300}-{tss + 300}"], text=True)
        for line in out.splitlines():
            fld = line.split("\t")
            flag = int(fld[1])
            pos = int(fld[3])
            L = cigar_len(fld[5])
            five, three = five_three(flag, pos, L)
            if strand == "+":
                prof5[five - tss] += 1
                prof3[three - tss] += 1
            else:
                prof5[tss - five] += 1
                prof3[tss - three] += 1

    print("\n========== (2) TSS 相对位置 profile（下游为正，多基因合并） ==========")
    print(f"{'relPos':>7s} {'5prime端':>8s} {'3prime端':>8s}")
    for r in range(-200, 301):
        if prof5.get(r, 0) or prof3.get(r, 0):
            print(f"{r:7d} {prof5.get(r, 0):8d} {prof3.get(r, 0):8d}")

    # 暂停区 +20~+100bp 内的 5' / 3' 端总计数，判断哪个端是聚合酶位点。
    pause5 = sum(prof5.get(r, 0) for r in range(20, 101))
    pause3 = sum(prof3.get(r, 0) for r in range(20, 101))
    peak5 = max(range(20, 101), key=lambda r: prof5.get(r, 0))
    peak3 = max(range(20, 101), key=lambda r: prof3.get(r, 0))
    print(f"\n暂停区(下游+20~+100bp): 5'端={pause5} (峰 rel={peak5})  3'端={pause3} (峰 rel={peak3})")
    if pause5 >= pause3:
        print("结论: 5' 端在暂停区富集 → signal_end = 5（read 5' 端 = 聚合酶位点）")
    else:
        print("结论: 3' 端在暂停区富集 → signal_end = 3（read 3' 端 = 聚合酶位点）")


if __name__ == "__main__":
    main()
