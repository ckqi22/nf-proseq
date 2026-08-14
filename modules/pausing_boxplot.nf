process pausing_boxplot {
    tag "pausing_boxplot"

    input:
    path pi_table          // PI table from pausing_index.R
    val  groups_yaml       // YAML text (group_name -> [samples])

    output:
    path "Pausing_Index_Boxplot.pdf", emit: pdf

    script:
    """
    source /workplace/hanguojun/mambaforge/bin/activate snakemake

    cat > pause_groups.yml << 'EOF'
    ${groups_yaml}
    EOF

    ${params.r} ${projectDir}/bin/pausing_boxplot.R \\
        --pi ${pi_table} \\
        --groups pause_groups.yml \\
        --output Pausing_Index_Boxplot.pdf
    """
}
