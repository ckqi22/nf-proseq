process DESEQ2 {
    tag "deseq2"

    input:
    path counts          // gene body count matrix (gene_id, Length, samples)
    val  config_yml      // differential config YAML (group + compared-groups)
    path annotation      // gene annotation (gene_id, gene_name, gene_biotype)

    output:
    path "deseq2_out/", emit: result

    script:
    """
    cat > deseq2_config.yml << 'EOF'
${config_yml}
EOF

    mkdir -p deseq2_out

    /workplace/shuixin/R/4.4.1/bin/Rscript /workplace/pipeline/WTSS/scripts/DESeq2_DE_analyses.R \\
        -f ${counts} \\
        -f2 ${counts} \\
        -a ${annotation} \\
        --config deseq2_config.yml \\
        --trans GeneBody \\
        -o deseq2_out
    """
}
