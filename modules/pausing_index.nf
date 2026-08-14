process pausing_index {
    tag "pausing_index"

    input:
    path tss_counts           // TSS count matrix (gene_id, length, <samples>...)
    path gene_body_counts     // gene body count matrix

    output:
    path "Pausing_Index.tsv", emit: pi

    script:
    def min_gene_length = params.tss.min_gene_length ?: 800
    """
    source /workplace/hanguojun/mambaforge/bin/activate snakemake

    ${params.r} ${projectDir}/bin/pausing_index.R \\
        --tss_counts ${tss_counts} \\
        --gb_counts ${gene_body_counts} \\
        --min_gene_length ${min_gene_length} \\
        --output Pausing_Index.tsv
    """
}
