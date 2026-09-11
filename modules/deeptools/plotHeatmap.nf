process PLOTHEATMAP {
    tag "${meta.sample}"

    conda '/workplace/hanguojun/mambaforge/envs/deeptools'

    input:
    tuple val(meta), path(matrix)

    output:
    tuple val(meta), path("${meta.sample}_metagene_heatmap.{pdf,tiff}"), emit: plot
    tuple val(meta), path("${meta.sample}_metagene_heatmap_matrix.gz"), emit: matrix
    tuple val(meta), path("${meta.sample}_metagene_heatmap_sorted_regions.bed"), emit: sorted_regions

    script:
    // 逐基因 TSS 热图（每行一个基因），按窗口均值升序，即弱信号在顶部。
    // 色阶：Reds（白→红）+ --zMin 0（矩阵值经 --scale 翻正后全 ≥0）。
    // （.gz，按排序后行序写出，与热图行序一致）。
    """
    plotHeatmap \\
        --matrixFile ${matrix} \\
        --sortRegions ascend \\
        --whatToShow "heatmap and colorbar" \\
        --colorMap Reds \\
        --zMin 0 \\
        --dpi 300 \\
        --heatmapHeight 20 \\
        --xAxisLabel "Distance from TSS (bp)" \\
        --outFileName ${meta.sample}_metagene_heatmap.pdf \\
        --outFileNameMatrix ${meta.sample}_metagene_heatmap_matrix.gz \\
        --outFileSortedRegions ${meta.sample}_metagene_heatmap_sorted_regions.bed

    convert \\
        -density 300 -quality 100 \\
        ${meta.sample}_metagene_heatmap.pdf \\
        ${meta.sample}_metagene_heatmap.tiff
    """
}
