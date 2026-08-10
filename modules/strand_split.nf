process strand_split {
    tag "${meta.sample}"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${meta.sample}.plus.bam"), path("${meta.sample}.minus.bam"), emit: strand_bams
    tuple val(meta), path("${meta.sample}.plus.bam.bai"), emit: plus_bai
    tuple val(meta), path("${meta.sample}.minus.bam.bai"), emit: minus_bai

    script:
    """
    source /workplace/hanguojun/mambaforge/bin/activate snakemake

    # Forward strand (-F 0x10: read NOT reverse, i.e. forward orientation)
    samtools view -@ 4 -b -F 0x10 ${bam} | \\
        samtools sort -@ 4 -o ${meta.sample}.plus.bam
    samtools index ${meta.sample}.plus.bam

    # Reverse strand (-f 0x10: read IS reverse orientation)
    samtools view -@ 4 -b -f 0x10 ${bam} | \\
        samtools sort -@ 4 -o ${meta.sample}.minus.bam
    samtools index ${meta.sample}.minus.bam
    """
}
