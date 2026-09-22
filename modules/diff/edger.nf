process EDGER {
    tag "edger"

    input:
    path counts          // gene body count matrix (gene_id, length, samples)
    val  config_yml      // differential config YAML (group + compared-groups)
    path annotation      // gene annotation (gene_id, gene_name, gene_biotype)

    output:
    path "edger_out/", emit: result

    script:
    """
    # The WTSS edgeR script computes FPKM from a 'Length' column (capital L);
    # our merged matrix uses 'length' (lowercase). Rename the header only.
    sed '1s/\\blength\\b/Length/' ${counts} > gene_body_counts.txt

    cat > edger_config.yml << 'EOF'
${config_yml}
EOF

    mkdir -p edger_out

    ${params.r} /workplace/pipeline/WTSS/scripts/edgeR_DE_analyses.R \\
        -f gene_body_counts.txt \\
        -f2 gene_body_counts.txt \\
        -a ${annotation} \\
        --config edger_config.yml \\
        --trans GeneBody \\
        -o edger_out
    """
}
