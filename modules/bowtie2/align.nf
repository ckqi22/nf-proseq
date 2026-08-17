process BOWTIE2_ALIGN {
    tag "${meta.sample}"

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), val(index)


    output:
    tuple val(meta), path("${meta.sample}.bam"), emit: bam
    tuple val(meta), path("${meta.sample}.bam.bai"), emit: bai
    tuple val(meta), path("${meta.sample}_summary_bowtie2.txt"), emit: alignRate


    script:
    def reads_args = ""
    def strand_args = ""

    if (meta.single_end) {
        reads_args = "-U ${reads}"
        // 单端用 --nofw / --norc
        if (params.strandedness == 'forward')       strand_args = "--norc"
        if (params.strandedness == 'reverse')       strand_args = "--nofw"
    } else {
        reads_args = "-1 ${reads[0]} -2 ${reads[1]}"
        // 双端用 --fr / --rf / --ff
        if (params.strandedness == 'forward')   strand_args = "--fr"
        if (params.strandedness == 'reverse')   strand_args = "--rf"
        if (params.strandedness == 'unstranded') strand_args = "--fr"
    }
    """
    source /workplace/hanguojun/mambaforge/bin/activate snakemake

    # NOTE: --no-discordant / --no-mixed 已禁用（原命令里这两个 flag 会被触发
    #   bowtie2 输出 SEQ/QUAL 长度不一致的 SAM 记录，导致
    #   `samtools view` 报 [E::sam_parse1] SEQ and QUAL are of different length。
    #   需要恢复时把下面两个 flag 加回 --no-unal 之后即可。）
    bowtie2 \\
        -x ${index} \\
        ${reads_args} \\
        --threads 10 \\
        ${strand_args} \\
        --very-sensitive \\
        --no-unal \\
        2> >(tee ${meta.sample}_summary_bowtie2.txt) | \\
    samtools view  -@ 2 -bS - | \\
    samtools sort  -@ 10 -o ${meta.sample}.bam && \\
    samtools index -@ 1 ${meta.sample}.bam
    """
}