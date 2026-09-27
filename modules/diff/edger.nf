process EDGER {
    tag "edger"

    input:
    path counts          // gene body count matrix (gene_id, Length, samples)
    val  config_yml      // differential config YAML (group + compared-groups)
    path annotation      // gene annotation (gene_id, gene_name, gene_biotype)

    output:
    path "edger_out/", emit: result

    script:
    """
    cat > edger_config.yml << 'EOF'
${config_yml}
EOF

    mkdir -p edger_out

    /workplace/shuixin/R/4.4.1/bin/Rscript /workplace/pipeline/WTSS/scripts/edgeR_DE_analyses.R \\
        -f ${counts} \\
        -f2 ${counts} \\
        -a ${annotation} \\
        --config edger_config.yml \\
        --trans GeneBody \\
        -o edger_out
    """
}
