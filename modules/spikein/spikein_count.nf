process SPIKEIN_COUNT {
    tag "${meta.sample}"

    // 从「主基因组 + spike 合并参考」比对出的 BAM 中，按 spike 染色体统计 read1 数。
    // 只数 read1（与下游信号口径一致）：SE 无 read2 即全量；PE 用 -f 64 取 first-in-pair。
    //   不用 idxstats（其列 3 按 read 段计，PE 会把 read1+read2 都算进去，分子偏大）；samtools view 无需 .bai。
    // spike_chroms 由 SPIKEIN_CONCAT（spikein_concat.nf）生成（spike fasta 头第一 token，逐行）。

    input:
    tuple val(meta), path(bam)
    path spike_chroms

    output:
    tuple val(meta), path("${meta.sample}.spike_count.txt"), emit: counts   // 文件内单行整数

    script:
    def read1_flag = meta.single_end ? '' : '-f 64 '
    """
    samtools view ${read1_flag}-F 4 ${bam} \\
        | awk 'NR==FNR{sp[\$1]=1; next} (\$3 in sp){s++} END{print s+0}' ${spike_chroms} - \\
        > ${meta.sample}.spike_count.txt
    """
}
