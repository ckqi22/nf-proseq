process BOWTIE2_ALIGN {
    tag "${meta.sample}"

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), val(index)


    output:
    tuple val(meta), path("${meta.sample}.bam"), emit: bam
    tuple val(meta), path("${meta.sample}.bam.bai"), emit: bai
    tuple val(meta), path("${meta.sample}_summary_bowtie2.txt"), emit: genomeRate
    tuple val(meta), path("z.${meta.sample}_bowtie2_samtools.log"), emit: alignment_log


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

    bowtie2 \\
        -x ${index} \\
        ${reads_args} \\
        --threads 10 \\
        ${strand_args} \\
        --very-sensitive \\
        --no-unal \\
        --no-discordant \\
        --no-mixed \\
        2> >(tee ${meta.sample}_summary_bowtie2.txt) | \\
    samtools view  --threads 10 -bS - | \\
    samtools sort  --threads 8 -o ${meta.sample}.bam && \\
    samtools index --threads 1 ${meta.sample}.bam > z.${meta.sample}_bowtie2_samtools.log 2>&1
    """
}