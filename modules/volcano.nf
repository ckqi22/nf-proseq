process volcano {
    tag "volcano"

    input:
    path de_results
    val config

    output:
    path "Volcano_Plot.pdf", emit: volcano_pdf

    script:
    def fc_cutoff   = params.diff?.fc_cutoff ?: params.diff_fc_cutoff ?: 1.5
    def pval_cutoff = params.diff?.pval_cutoff ?: params.diff_pval_cutoff ?: 0.05

    """
    source /workplace/hanguojun/mambaforge/bin/activate snakemake

    ${params.r} ${projectDir}/bin/volcano.R \\
        --de_results ${de_results} \\
        --fc_cutoff ${fc_cutoff} \\
        --pval_cutoff ${pval_cutoff} \\
        --output Volcano_Plot.pdf
    """
}
