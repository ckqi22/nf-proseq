process geo_prep {
    tag "geo_prep"

    input:
    path samplesheet
    path counts_matrix
    val config

    output:
    path "GEO_submission", emit: geo_output

    script:
    """
    source /workplace/hanguojun/mambaforge/bin/activate snakemake

    mkdir -p GEO_submission

    # Generate GEO submission templates
    # Processed data files, metadata spreadsheet, and README

    python3 ${projectDir}/bin/geo_submission_prep.py \\
        --samplesheet ${samplesheet} \\
        --counts_matrix ${counts_matrix} \\
        --output_dir GEO_submission
    """
}
