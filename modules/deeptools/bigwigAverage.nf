process BIGWIGAVERAGE {
    tag "${meta.group}"

    input:
    tuple val(meta), path(plus_bws), path(minus_bws)
    path chrom_sizes

    output:
    tuple val(meta), path("${meta.group}${meta.plus_suffix}.bigWig"), path("${meta.group}${meta.minus_suffix}.bigWig"), emit: avg_bigwig

    script:
    // 组内样本分析 bigWig 平均（重复合并）：
    // 输入应该已归一化（cpm/spike）→ 均值即组内平均。
    // minus 链是负值：bigWigMerge 默认 -threshold=0.0 会丢弃 ≤0 的区间，故必须 -threshold=-1e30（足够负）保留负值。
    // 单样本的组不调工具，直接复制（平均值 = 自身）。
    // 排序：bigWigMerge 按 bigWig 内部 bbi 索引顺序输出（实测自然序 chr1..chr22,chrM,chrX,chrY），
    // 不是 bedGraphToBigWig v2.8 要的 C 字典序（chr1,chr10,..,chr19,chr2,..），
    // 故平均后的 bedGraph 须 LC_COLLATE=C sort 重排再写回；chrom.sizes 顺序无关（工具只查大小，不查顺序）。
    def plus_out  = "${meta.group}${meta.plus_suffix}.bigWig"
    def minus_out = "${meta.group}${meta.minus_suffix}.bigWig"
    def plus_n    = plus_bws.size()
    def minus_n   = minus_bws.size()
    def plus_cmd = plus_n > 1
        ? "bigWigMerge ${plus_bws.join(' ')}  plus.merged.bedGraph && " +
          "awk -v n=${plus_n} 'BEGIN{OFS=\"\\t\"}{\$4=\$4/n; if(\$4!=0) print}' plus.merged.bedGraph | LC_COLLATE=C sort -k1,1 -k2,2n > plus.avg.bedGraph && " +
          "bedGraphToBigWig plus.avg.bedGraph chrom.sizes ${plus_out}"
        : "cp ${plus_bws[0]} ${plus_out}"
    def minus_cmd = minus_n > 1
        ? "bigWigMerge -threshold=-1e30 ${minus_bws.join(' ')} minus.merged.bedGraph && " +
          "awk -v n=${minus_n} 'BEGIN{OFS=\"\\t\"}{\$4=\$4/n; if(\$4!=0) print}' minus.merged.bedGraph | LC_COLLATE=C sort -k1,1 -k2,2n > minus.avg.bedGraph && " +
          "bedGraphToBigWig minus.avg.bedGraph chrom.sizes ${minus_out}"
        : "cp ${minus_bws[0]} ${minus_out}"
    """
    cut -f1,2 ${chrom_sizes} > chrom.sizes
    ${plus_cmd}
    ${minus_cmd}
    """
}


// ---- 原 bigwigAverage 实现（备选/对照，已停用）----
// def plus_cmd = plus_bws.size() > 1
//     ? "bigwigAverage -b ${plus_bws.join(' ')}  --skipNAs -bs 1 -p 40 -o ${plus_out}"
//     : "cp ${plus_bws[0]} ${plus_out}"
// def minus_cmd = minus_bws.size() > 1
//     ? "bigwigAverage -b ${minus_bws.join(' ')} --skipNAs -bs 1 -p 40 -o ${minus_out}"
//     : "cp ${minus_bws[0]} ${minus_out}"
// 说明：bigwigAverage 按 -bs 分箱后取 bin 内均值，-bs 1 逐碱基分箱遍历全基因组，
//   速度极慢
