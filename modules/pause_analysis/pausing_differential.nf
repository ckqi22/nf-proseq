process pausing_differential {
    tag "pausing_differential"

    container "bio-base:1.0.0"

    input:
    path pi_table          // PI table from pausing_index.R
    val  groups_yaml       // YAML text (group_name -> [samples])
    val  comparisons_csv   // CSV text (group1,group2)

    output:
    path "Pausing_Differential.tsv", emit: diff

    script:
    """
    source /home/ck/miniconda3/bin/activate renv

    cat > pause_groups.yml << 'EOF'
    ${groups_yaml}
    EOF

    cat > pause_comparisons.csv << 'EOF'
    ${comparisons_csv}
    EOF

    Rscript ${projectDir}/bin/pausing_differential.R \\
        --pi ${pi_table} \\
        --groups pause_groups.yml \\
        --comparisons pause_comparisons.csv \\
        --output Pausing_Differential.tsv
    """
}
