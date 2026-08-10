process pausing_index {
    tag "pausing_index"

    input:
    path tss_counts           // Combined TSS count matrix (all samples)
    path gene_body_counts     // Combined gene body count matrix (all samples)
    path groups_yaml          // Optional: YAML file mapping group names to sample lists
    path comparisons_csv      // Optional: CSV with columns: group1, group2

    output:
    path "Pausing_Index_All.xlsx",         emit: pi_all
    path "Pausing_Index_Boxplot.pdf",      emit: pi_boxplot
    path "Pausing_Differential.xlsx",      emit: pi_diff

    script:
    def min_gene_length  = params.min_gene_length ?: 800
    def groups_arg       = (groups_yaml && groups_yaml.name != 'null') ? "--groups ${groups_yaml}" : ""
    def comp_arg         = (comparisons_csv && comparisons_csv.name != 'null') ? "--comparisons ${comparisons_csv}" : ""

    """
    source /workplace/hanguojun/mambaforge/bin/activate snakemake

    ${params.r} ${projectDir}/bin/pausing_index.R \\
        --tss_counts ${tss_counts} \\
        --gb_counts ${gene_body_counts} \\
        --min_gene_length ${min_gene_length} \\
        ${groups_arg} \\
        ${comp_arg} \\
        --output_dir ./

    test -f Pausing_Index_All.xlsx || { echo "ERROR: Pausing_Index_All.xlsx not generated"; exit 1; }
    """
}
