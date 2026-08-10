process edger_de {
    tag "edger_de"

    input:
    path counts_matrix
    path groups_yaml         // YAML: group_name → [sample1, sample2, ...]
    path comparisons_csv     // CSV: group1, group2
    val config

    output:
    path "Diff_genes.xlsx",              emit: diff_genes
    path "All_Comparisons_genes.xlsx",   emit: all_comparisons

    script:
    def fc_cutoff   = params.diff?.fc_cutoff ?: params.diff_fc_cutoff ?: 1.5
    def pval_cutoff = params.diff?.pval_cutoff ?: params.diff_pval_cutoff ?: 0.05

    """
    source /workplace/hanguojun/mambaforge/bin/activate snakemake

    ${params.r} ${projectDir}/bin/edger_de.R \\
        --counts ${counts_matrix} \\
        --groups ${groups_yaml} \\
        --comparisons ${comparisons_csv} \\
        --fc_cutoff ${fc_cutoff} \\
        --pval_cutoff ${pval_cutoff} \\
        --output_dir ./
    """
}
