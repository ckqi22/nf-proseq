process PLOTHEATMAP {
    tag "${meta.sample}"

    input:
    tuple val(meta), path(matrix)

    output:
    tuple val(meta), path("${meta.sample}_metagene_heatmap.pdf"), emit: plot

    script:
    """
    source /workplace/hanguojun/mambaforge/bin/activate deeptools

    plotProfile \\
        --matrixFile ${matrix} \\
        --outFileName ${meta.sample}_TSS_meta.pdf \\
        --sortRegions ascend \\
        --whatToShow "heatmap and colorbar"
    """
}
