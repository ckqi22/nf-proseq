process DESEQ2 {
    tag "deseq2"

    // Runs on the host (no container) so it can reach ${params.r} and the fixed
    // WTSS script path below.
    input:
    path counts          // gene body count matrix (gene_id, length, samples)
    val  config_yml      // differential config YAML (group + compared-groups)
    path annotation      // gene annotation (gene_id, gene_name, gene_biotype)

    output:
    path "deseq2_out/*", emit: results

    script:
    """
    # The WTSS DESeq2 script computes FPKM from a 'Length' column (capital L);
    # our merged matrix uses 'length' (lowercase). Rename the header only.
    sed '1s/\\blength\\b/Length/' ${counts} > gene_body_counts.txt

    cat > deseq2_config.yml << 'EOF'
${config_yml}
EOF

    mkdir -p deseq2_out

    ${params.r} /workplace/pipeline/WTSS/scripts/DESeq2_DE_analyses.R \\
        -f gene_body_counts.txt \\
        -f2 gene_body_counts.txt \\
        -a ${annotation} \\
        --config deseq2_config.yml \\
        --trans GeneBody \\
        -o deseq2_out
    """
}
