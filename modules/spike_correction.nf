process spike_correction {
    tag "spike_correction"

    input:
    path spike_counts
    path gene_body_counts
    tuple val(meta)

    output:
    path "spike_correction_factors.txt", emit: spike_factors

    script:
    """
    source /workplace/hanguojun/mambaforge/bin/activate snakemake

    # Compute spike-in correction factors
    # spike_counts: featureCounts output aligned to spike genome (dm6)
    # gene_body_counts: endogenous gene body counts
    # Output: sigma_j (scale factor per sample) and epsilon_j (normalization factor)

    ${params.r} ${projectDir}/bin/spike_correction.R \\
        --spike ${spike_counts} \\
        --endogenous ${gene_body_counts} \\
        --output spike_correction_factors.txt
    """
}
