process scatter {
    tag "scatter"

    input:
    path rpkm_matrix
    val comparisons
    val config

    output:
    path "Scatter_Plot.pdf", emit: scatter_pdf

    script:
    """
    source /workplace/hanguojun/mambaforge/bin/activate snakemake

    # Scatter plot of sample correlations using ggplot2

    ${params.r} ${projectDir}/bin/scatter.R \\
        --rpkm ${rpkm_matrix} \\
        --comparisons "${comparisons}" \\
        --output Scatter_Plot.pdf
    """
}
