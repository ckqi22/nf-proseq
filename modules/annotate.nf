process ANNOTATE {
    tag "${name}"

    container "nf-proseq:1.0.0"

    input:
    path raw             // raw count matrix (gene_id, length, samples)
    path normalized      // normalize.R output (gene_id, <s>.<Suffix>)
    path annotation      // local gene annotation table (first col gene_id)
    val name             // profile name (genebody | promoter | pol2)

    output:
    path "${name}.annotated.txt", emit: annotated

    script:
    """
    Rscript ${projectDir}/bin/annotate.R \\
        --raw ${raw} \\
        --normalized ${normalized} \\
        --annotation ${annotation} \\
        --output ${name}.annotated.txt
    """
}
