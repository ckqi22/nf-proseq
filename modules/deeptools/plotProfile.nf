process PLOTPROFILE {
    tag "${meta.sample}"

    conda '/workplace/hanguojun/mambaforge/envs/deeptools'

    input:
    tuple val(meta), path(matrix)

    output:
    tuple val(meta), path("${meta.sample}_metagene_profile_${meta.sig}.{pdf,tiff}"), emit: profile
    tuple val(meta), path("${meta.sample}_metagene_profile_${meta.sig}_matrix.tsv"), emit: matrix

    script:
    // --perGroup：叠加图必加（metagene.nf / metagene_group.nf 在 meta 写 per_group:true）。
    // 不加时 plotProfile 默认 numplots=样本数 → 每个样本各占一个 panel；加了之后
    // numplots=region group 数（rbind 已合并为单 group=1）、numlines=样本数 → 一 panel 多线。
    def pg = meta.per_group ? '--perGroup' : ''
    // y 轴标签随信号描述符 sig（single_cpm / single_spike / full_cpm / full_spike）
    def ylabel = [
        single_cpm:   "5' end CPM",
        single_spike: "5' end RPM (spike)",
        full_cpm:     "full-read CPM",
        full_spike:   "full-read RPM (spike)"
    ].get(meta.sig, "CPM")
    """
    plotProfile \\
        --matrixFile ${matrix} \\
        --plotType lines \\
        ${pg} \\
        --dpi 300 \\
        --yAxisLabel "${ylabel}" \\
        --outFileName ${meta.sample}_metagene_profile_${meta.sig}.pdf \\
        --outFileNameData ${meta.sample}_metagene_profile_${meta.sig}_matrix.tsv

    convert \\
        -density 300 -quality 100 \\
        ${meta.sample}_metagene_profile_${meta.sig}.pdf \\
        ${meta.sample}_metagene_profile_${meta.sig}.tiff
    """
}
