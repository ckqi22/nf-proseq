process GTF2SAF {

    container "bio-base:1.0.0"

    input:
    path gtf

    output:
    tuple val("genebody"), path('genebody_union.saf'), emit: genebody_union_saf // 所有 transcript genebody union → quantification(featureCounts)

    script:
    """
    source /home/ck/miniconda3/bin/activate renv
    Rscript ${projectDir}/bin/gtf2saf.R \\
        --gtf ${gtf} \\
        --genebody_offset ${params.tss.genebody_offset} \\
        --outdir ./
    """
}
