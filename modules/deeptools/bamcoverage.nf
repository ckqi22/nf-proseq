process BAMCOVERAGE {
    tag "${meta.sample}"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${meta.sample}.bw"), emit: bigwig

    script:
    """
    source /workplace/hanguojun/mambaforge/bin/activate deeptools

    bamCoverage \\
        --bam ${bam} \\
        --outFileName ${meta.sample}.bw \\
        --outFileFormat bigwig \\
        --binSize 1 \\
        --Offset 1 \\
        --numberOfProcessors 10 \\
        --normalizeUsing RPKM
    """
}
