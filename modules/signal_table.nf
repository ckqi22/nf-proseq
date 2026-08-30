process SIGNAL_TABLE {
    tag "signal_table"

    container "bio-base:1.0.0"

    input:
    path bedgraph_files      // flat list of *_plus.bedgraph / *_minus.bedgraph（staged，供 R 按 basename 读取）
    val manifest              // "sample \t plus_bg \t minus_bg"（样本名来自 meta.sample）
    val groups_yaml           // YAML 字符串：group -> [samples]
    path gene_bed             // 标准 BED6
    path gtf                  // 参考 GTF（取代表 transcript）
    path annotation           // gene 注释表（首列 gene_id）

    output:
    path "pol2_signal_table.tsv", emit: signal_table

    script:
    """
    source /home/ck/miniconda3/bin/activate renv

    cat > manifest.tsv << 'EOF'
${manifest}
EOF

    cat > groups.yml << 'EOF'
${groups_yaml}
EOF

    Rscript ${projectDir}/bin/signal_table.R \\
        --manifest manifest.tsv \\
        --groups groups.yml \\
        --gene_bed ${gene_bed} \\
        --gtf ${gtf} \\
        --annotation ${annotation} \\
        --output pol2_signal_table.tsv
    """
}
