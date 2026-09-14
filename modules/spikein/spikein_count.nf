process SPIKEIN_COUNT {
    tag "${meta.sample}"

    // 从「主基因组 + spike 合并参考」比对出的 BAM 中，按 spike 染色体统计 read 数。
    // samtools idxstats 依赖 .bai，故输入连同 .bai 一起 stage（文件名须为 ${bam}.bai）。
    // spike_chroms 由 SPIKEIN_CONCAT（spikein_concat.nf）生成（spike fasta 头第一 token，逐行）。

    input:
    tuple val(meta), path(bam), path(bai)
    path spike_chroms

    output:
    tuple val(meta), path("${meta.sample}.spike_count.txt"), emit: counts   // 文件内单行整数

    script:
    """
    samtools idxstats ${bam} \\
        | awk 'NR==FNR{sp[\$1]=1; next} (\$1 in sp){s+=\$3} END{print s+0}' ${spike_chroms} - \\
        > ${meta.sample}.spike_count.txt
    """
}
