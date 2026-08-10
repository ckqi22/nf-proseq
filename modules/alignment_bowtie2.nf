process alignment_bowtie2 {
    tag "${meta.sample}"

    input:
    tuple val(meta), path(reads)
    val config


    output:
    tuple val(meta), path("${meta.sample}.bam"), emit: bam
    tuple val(meta), path("${meta.sample}.bam.bai"), emit: bai
    tuple val(meta), path("${meta.sample}_summary_bowtie2.txt"), emit: genomeRate
    tuple val(meta), path("z.${meta.sample}_hisat2_samtools.log"), emit: alignment_log


    script:
    """
    source /workplace/hanguojun/mambaforge/bin/activate snakemake
    
    bowtie2 \\
    -x ${config.bowtie2_index} \\
    -1 ${reads[0]} -2 ${reads[1]} \\
    -p 10 \\
    -q \\
    --rf \\
    --very-sensitive \\
    --no-unal \\
    --no-discordant \\
    --no-mixed \\
    2>${meta.sample}_summary_bowtie2.txt | \\
    samtools view -@ 10 -bS - | \\
    samtools sort -@ 8 -o ${meta.sample}.bam && \\
    samtools index ${meta.sample}.bam > z.${meta.sample}_hisat2_samtools.log 2>&1
    """
}
