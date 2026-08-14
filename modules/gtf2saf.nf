process GTF2SAF {

    container "bio-base:1.0.0"

    input:
    path gtf

    output:
    path '*.saf', emit: saf

    script:
    """
    source /home/ck/miniconda3/bin/activate renv
    Rscript ${projectDir}/bin/gtf2saf.R \\
        --gtf ${gtf} \\
        --tss_upstream ${params.tss.upstream} \\
        --tss_downstream ${params.tss.downstream} \\
        --genebody_offset ${params.tss.genebody_offset} \\
        --outdir ./
    """
}