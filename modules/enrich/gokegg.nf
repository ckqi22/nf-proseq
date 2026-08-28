process GOKEGG {

    input:
    path diff_dir

    output:
    path "gokegg_result/*", emit: result

    script:
    """
    mkdir gokegg_result

    ${params.r} /workplace/pipeline/code/enrichment.R \\
        --species ${params.species} \\
        --input_dir ${diff_dir} \\
        --output_dir ./gokegg_result \\
        --methods NORMAL \\
        --Pgo 1 --Pkegg 1 --PKEGGgsea 1 --PGOgsea 1
    """
}