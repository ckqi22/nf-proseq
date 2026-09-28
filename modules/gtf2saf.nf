process GTF2SAF {

    container "nf-proseq:1.0.0"

    input:
    path gtf

    output:
    tuple val("genebody"), path('genebody_union.saf'), emit: genebody_union_saf // 所有 transcript genebody union → quantification(featureCounts)

    script:
    """
    Rscript ${projectDir}/bin/gtf2saf.R \\
        --gtf ${gtf} \\
        --genebody_offset ${params.tss.genebody_offset} \\
        --outdir ./
    """
}
