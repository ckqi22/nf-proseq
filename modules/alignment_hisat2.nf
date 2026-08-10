process alignment {
    tag "${meta.sample}"

    input:
    tuple val(meta), path(reads)
    val config


    output:
    tuple val(meta), path("${meta.sample}.bam"), emit: bam
    tuple val(meta), path("${meta.sample}.bam.bai"), emit: bai
    tuple val(meta), path("${meta.sample}_summary_hisat2.txt"), emit: genomeRate
    tuple val(meta), path("z.${meta.sample}_hisat2_samtools.log"), emit: alignment_log


    script:
    """
    source /workplace/hanguojun/mambaforge/bin/activate snakemake
    
    ${params.hisat2} \\
    -p 10 \\
    -q \\
    --rna-strandness RF \\
    -x ${config.hisat2_index} \\
    --summary-file ${meta.sample}_summary_hisat2.txt \\
    -1 ${reads[0]} -2 ${reads[1]} | \\
    samtools view -bS -F 4 -F 8 -F 256 - | \\
    samtools sort -@ 8 -o ${meta.sample}.bam && samtools index ${meta.sample}.bam > z.${meta.sample}_hisat2_samtools.log 2>&1
    """
}