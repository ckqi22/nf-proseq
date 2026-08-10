process bigwig_gen {
    tag "${meta.sample}"

    input:
    tuple val(meta), path(plus_bam), path(minus_bam)

    output:
    tuple val(meta), path("${meta.sample}.plus.bw"), path("${meta.sample}.minus.bw"), emit: bigwigs

    script:
    def effective_genome_size = params.effective_genome_size ?: 2913022398

    """
    source /workplace/hanguojun/mambaforge/bin/activate snakemake

    # Generate strand-specific BigWig files using deepTools bamCoverage
    # RPGC normalization for comparable tracks

    bamCoverage \\
        -b ${plus_bam} \\
        -o ${meta.sample}.plus.bw \\
        --normalizeUsing RPGC \\
        --effectiveGenomeSize ${effective_genome_size} \\
        --binSize 10 \\
        --smoothLength 30 \\
        -p 4

    bamCoverage \\
        -b ${minus_bam} \\
        -o ${meta.sample}.minus.bw \\
        --normalizeUsing RPGC \\
        --effectiveGenomeSize ${effective_genome_size} \\
        --binSize 10 \\
        --smoothLength 30 \\
        -p 4
    """
}
