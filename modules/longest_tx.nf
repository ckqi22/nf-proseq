process LONGEST_TX {

    container "bio-base:1.0.0"

    input:
    path gtf

    output:
    path 'longest_tx.gtf', emit: representative_gtf

    script:
    """
    source /home/ck/miniconda3/bin/activate renv
    Rscript ${projectDir}/bin/longest_tx.R \\
        --gtf ${gtf} \\
        --gene-type protein_coding \\
        --outdir ./
    """
}
