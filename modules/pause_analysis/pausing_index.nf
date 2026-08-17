process pausing_index {
    tag "pausing_index"

    container "bio-base:1.0.0"

    input:
    path tss_counts           // TSS count matrix (gene_id, length, <samples>...)
    path genebody_counts     // gene body count matrix

    output:
    path "Pausing_Index.tsv", emit: pi

    script:
    def min_gene_length = params.tss.min_gene_length ?: 800
    """
    source /home/ck/miniconda3/bin/activate renv

    Rscript ${projectDir}/bin/pausing_index.R \\
        --tss_counts ${tss_counts} \\
        --gb_counts ${genebody_counts} \\
        --min_gene_length ${min_gene_length} \\
        --output Pausing_Index.tsv
    """
}
