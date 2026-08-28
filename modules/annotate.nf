process ANNOTATE {
    tag "${name}"

    container "bio-base:1.0.0"

    input:
    path raw             // raw count matrix (gene_id, length, samples)
    path normalized      // normalize.R output (gene_id, <s>.<Suffix>)
    path annotation      // local gene annotation table (first col gene_id)
    val name             // profile name (genebody | promoter | pol2)

    output:
    path "${name}.annotated.txt", emit: annotated

    script:
    """
    source /home/ck/miniconda3/bin/activate renv

    Rscript ${projectDir}/bin/annotate.R \\
        --raw ${raw} \\
        --normalized ${normalized} \\
        --annotation ${annotation} \\
        --output ${name}.annotated.txt
    """
}
