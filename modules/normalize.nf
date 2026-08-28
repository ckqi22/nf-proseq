process NORMALIZE {
    tag "${name}"

    container "bio-base:1.0.0"

    input:
    path matrix          // count matrix (gene_id, length, samples)
    val methods          // comma-separated: cpm,fpkm[,rpkm]
    val name             // profile name (genebody | promoter | pol2)

    output:
    path "${name}.normalized.txt", emit: normalized

    script:
    """
    source /home/ck/miniconda3/bin/activate renv

    Rscript ${projectDir}/bin/normalize.R \\
        --input ${matrix} \\
        --methods ${methods} \\
        --output ${name}.normalized.txt
    """
}
