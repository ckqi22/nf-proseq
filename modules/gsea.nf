process gsea {
    tag "gsea"

    input:
    path counts_matrix
    path groups_yaml       // YAML: group_name → [sample1, ...]
    path comparisons_csv   // CSV: group1, group2
    val config

    output:
    path "GSEA_results", emit: gsea_results

    script:
    """
    source /workplace/hanguojun/mambaforge/bin/activate snakemake

    mkdir -p GSEA_results

    ${params.r} ${projectDir}/bin/gsea.R \\
        --counts ${counts_matrix} \\
        --groups ${groups_yaml} \\
        --comparisons ${comparisons_csv} \\
        --output_dir GSEA_results
    """
}
