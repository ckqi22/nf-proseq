process COMPUTEMATRIX {
    tag "${meta.sample}"

    input:
    tuple val(meta), path(bw)
    path  gene_bed

    output:
    tuple val(meta), path("${meta.sample}_TSS_matrix.gz"), emit: matrix

    script:
    """
    source /workplace/hanguojun/mambaforge/bin/activate deeptools

    computeMatrix reference-point \\
        -R ${gene_bed} \\
        -S ${bw} \\
        --missingDataAsZero \\
        --referencePoint TSS \\
        --upstream 50 \\
        --downstream 150 \\
        -o ${meta.sample}_TSS_matrix.gz \\
        -p 20
    """
}

// ${meta.sample}_TSS_matrix.gz格式说明
// E.g: chr1    3889124 3900293 ENSG00000198912 .   -   1.000000    nan ...
// 每行对应一个基因组区域，第7列及之后各列对应一个计算窗口（bin）
// 计算窗口数由--upstream、--downstream 和 --binSize (默认10bp)共同决定
// 如TSS + 上下游各1000bp，则应有 1000/10 * 2 = 200 个计算窗口，计算窗口的值默认为mean
// 若有多个样本，如两个样本，则共有 6 + 200 * 2 = 406 列