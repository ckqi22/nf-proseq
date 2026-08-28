process GSEA {

    input:
    path diff_dir

    output:
    path "gsea_result/*", emit: result

    script:
    """
    mkdir gsea_result

    ${params.r} /workplace/pipeline/code/enrichment.R \\
        --species ${params.species} \\
        --input_dir ${diff_dir} \\
        --output_dir ./gsea_result \\
        --methods GSEA_KEGG \\
        --Pgo 1 --Pkegg 1 --PKEGGgsea 1 --PGOgsea 1
    """
}