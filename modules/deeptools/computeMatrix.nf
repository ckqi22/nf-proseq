process COMPUTEMATRIX {
    tag "${meta.sample}"

    input:
    tuple val(meta), path(plus_bw), path(minus_bw)
    path  gene_bed

    output:
    tuple val(meta), path("${meta.sample}_TSS_matrix.gz"), emit: matrix

    script:
    """
    source /workplace/hanguojun/mambaforge/bin/activate deeptools

    computeMatrix reference-point \\
        -R ${gene_bed} \\
        -S ${plus_bw} ${minus_bw} \\
        --referencePoint TSS \\
        --upstream ${params.metagene.tss_window} \\
        --downstream ${params.metagene.tss_window} \\
        -o ${meta.sample}_TSS_matrix.gz \\
        -p 20
    """
}
