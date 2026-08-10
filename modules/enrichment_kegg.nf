process enrichment_kegg {
    tag "kegg_enrichment"

    input:
    path de_results
    val config

    output:
    path "KEGG_results", emit: kegg_results

    script:
    """
    source /workplace/hanguojun/mambaforge/bin/activate snakemake

    mkdir -p KEGG_results

    # KEGG pathway enrichment analysis using clusterProfiler

    ${params.r} ${projectDir}/bin/enrichment_kegg.R \\
        --de_results ${de_results} \\
        --output_dir KEGG_results
    """
}
