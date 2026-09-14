process PLOTPROFILE {
    tag "${meta.sample}"

    conda '/workplace/hanguojun/mambaforge/envs/deeptools'

    input:
    tuple val(meta), path(matrix)

    output:
    tuple val(meta), path("${meta.sample}_metagene_profile.{pdf,tiff}"), emit: profile
    tuple val(meta), path("${meta.sample}_metagene_profile_matrix.tsv"), emit: matrix

    script:
    // plotType：单样本默认 se（均值±SE 带）；叠加图由 metagene.nf / metagene_group.nf 在
    // meta 写 plot_type:'lines' 得干净折线（一图多线，图例取自矩阵内 --samplesLabel）。
    // --perGroup：叠加图必加（metagene.nf / metagene_group.nf 在 meta 写 per_group:true）。
    // 不加时 plotProfile 默认 numplots=样本数 → 每个样本各占一个 panel；加了之后
    // numplots=region group 数（rbind 已合并为单 group=1）、numlines=样本数 → 一 panel 多线。
    def pt = meta.plot_type ?: 'se'
    def pg = meta.per_group ? '--perGroup' : ''
    """
    plotProfile \\
        --matrixFile ${matrix} \\
        --plotType ${pt} \\
        ${pg} \\
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
