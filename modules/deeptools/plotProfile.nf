process PLOTPROFILE {
    tag "${meta.sample}"

    conda '/workplace/hanguojun/mambaforge/envs/deeptools'

    input:
    tuple val(meta), path(matrix)

    output:
    tuple val(meta), path("${meta.sample}_metagene_profile.{pdf,tiff}"), emit: profile
    tuple val(meta), path("${meta.sample}_metagene_profile_matrix.tsv"), emit: matrix

    script:
    """
    plotProfile \\
        --matrixFile ${matrix} \\
        --plotType se \\
        --dpi 300 \\
        --yAxisLabel "5' end CPM" \\
        --outFileName ${meta.sample}_metagene_profile.pdf \\
        --outFileNameData ${meta.sample}_metagene_profile_matrix.tsv

    convert \\
        -density 300 -quality 100 \\
        ${meta.sample}_metagene_profile.pdf \\
        ${meta.sample}_metagene_profile.tiff
    """
}
