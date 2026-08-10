process heatmap {
    tag "heatmap"

    input:
    path de_results
    path rpkm_matrix
    val config

    output:
    path "Heatmap.pdf", emit: heatmap_pdf

    script:
    """
    source /workplace/hanguojun/mambaforge/bin/activate snakemake

    # Heatmap of differentially expressed genes using ComplexHeatmap

    ${params.r} ${projectDir}/bin/heatmap.R \\
        --de_results ${de_results} \\
        --rpkm ${rpkm_matrix} \\
        --output Heatmap.pdf
    """
}
