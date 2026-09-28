process pausing_index {
    tag "pausing_index"

    container "nf-proseq:1.0.0"

    input:
    path promoter_count     // promoter count matrix (gene_id, length, <samples>...)
    path genebody_count     // genebody count matrix

    output:
    path "Pausing_Index.tsv", emit: pi

    script:
    def min_genebody_length = params.tss.min_genebody_length ?: 800
    def pseudocount     = params.tss.pi_pseudocount ?: 1e-3
    """
    Rscript ${projectDir}/bin/pausing_index.R \\
        --promoter_count ${promoter_count} \\
        --genebody_count ${genebody_count} \\
        --min_genebody_len ${min_genebody_length} \\
        --pseudocount ${pseudocount} \\
        --output Pausing_Index.tsv
    """
}
