process COMPUTEMATRIX {
    tag "${meta.sample}"

    input:
    tuple val(meta), path(bigwigs)        // all samples' bigWigs (space-separated for -S)
    val gtf             // GTF path (TSS reference points)

    output:
    tuple val(meta), path("${meta.sample}_TSS_matrix.gz"), emit: matrix

    script:
    """
    source /workplace/hanguojun/mambaforge/bin/activate deeptools

    computeMatrix reference-point \\
        -R ${gtf} \\
        -S ${bigwigs} \\
        --referencePoint TSS \\
        --upstream ${params.metagene.tss_window} \\
        --downstream ${params.metagene.tss_window} \\
        -o ${meta.sample}_TSS_matrix.gz \\
        -p 15
    """
}
