process BIGWIGAVERAGE {
    tag "${meta.group}"

    conda '/workplace/hanguojun/mambaforge/envs/deeptools'

    input:
    tuple val(meta), path(plus_bws), path(minus_bws)

    output:
    tuple val(meta), path("${meta.group}_plus_cpm.bigWig"), path("${meta.group}_minus_cpm.bigWig"), emit: avg_bigwig

    script:
    // 组内样本的 5'-端 CPM bigWig 平均（重复合并）：
    // bigwigAverage 把基因组按 bin 切块后取各文件均值；-bs 10 按 10bp bin 平均。
    // 输入已 CPM 归一化 → 均值即组内平均 CPM。
    // 非覆盖区默认按 0 处理（与下游 computeMatrix --missingDataAsZero 一致）。
    // 单样本的组不调工具，直接复制（平均值 = 自身）。
    def plus_cmd = plus_bws.size() > 1
        ? "bigwigAverage -b ${plus_bws.join(' ')}  -bs 10 -p 20 -o ${meta.group}_plus_cpm.bigWig"
        : "cp ${plus_bws[0]} ${meta.group}_plus_cpm.bigWig"
    def minus_cmd = minus_bws.size() > 1
        ? "bigwigAverage -b ${minus_bws.join(' ')} -bs 10 -p 20 -o ${meta.group}_minus_cpm.bigWig"
        : "cp ${minus_bws[0]} ${meta.group}_minus_cpm.bigWig"
    """
    ${plus_cmd}
    ${minus_cmd}
    """
}
