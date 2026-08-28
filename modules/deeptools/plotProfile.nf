process PLOTPROFILE {
    tag "${meta.sample}"

    input:
    tuple val(meta), path(matrix)

    output:
    tuple val(meta), path("${meta.sample}_TSS_meta.pdf"), emit: profile

    script:
    """
    source /workplace/hanguojun/mambaforge/bin/activate deeptools

    plotProfile \\
        --matrixFile ${matrix} \\
        --outFileName ${meta.sample}_TSS_meta.pdf \\
        --perGroup \\
        --plotType lines \\
        --dpi 300
    """
}
