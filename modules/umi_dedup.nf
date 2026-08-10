process umi_dedup {
    tag "${meta.sample}"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${meta.sample}.dedup.bam"), emit: dedup_bam

    script:
    """
    source /workplace/hanguojun/mambaforge/bin/activate snakemake

    # UMI deduplication using umi_tools
    # The UMI was extracted from reads and prepended to the read name by fastp (--umi_prefix=UMI)
    umi_tools dedup \\
        --stdin=${bam} \\
        --stdout=${meta.sample}.dedup.bam \\
        --extract-umi-method=read_id \\
        --method=unique \\
        --log=${meta.sample}.umi_dedup.log

    samtools index ${meta.sample}.dedup.bam
    """
}
