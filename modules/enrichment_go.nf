process enrichment_go {
    tag "go_enrichment"

    input:
    path de_results
    val config

    output:
    path "GO_results", emit: go_results

    script:
    """
    source /workplace/hanguojun/mambaforge/bin/activate snakemake

    mkdir -p GO_results

    # GO enrichment analysis using clusterProfiler
    # BP (Biological Process), MF (Molecular Function), CC (Cellular Component)
    # For both up-regulated and down-regulated genes separately

    ${params.r} ${projectDir}/bin/enrichment_go.R \\
        --de_results ${de_results} \\
        --output_dir GO_results
    """
}
