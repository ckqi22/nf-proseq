process pca_plot {
    tag "pca"

    input:
    path counts_matrix
    val config

    output:
    path "PCA_Plot.pdf", emit: pca_pdf

    script:
    """
    source /workplace/hanguojun/mambaforge/bin/activate snakemake

    # PCA plot using R

    ${params.r} ${projectDir}/bin/pca_plot.R \\
        --counts ${counts_matrix} \\
        --output PCA_Plot.pdf
    """
}
