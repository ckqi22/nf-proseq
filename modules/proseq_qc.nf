process proseq_qc {
    tag "${meta.sample}"

    input:
    tuple val(meta), path(bam), path(alignment_log)

    output:
    path "${meta.sample}.PROseq_QC_report.txt", emit: qc_report

    script:
    """
    source /workplace/hanguojun/mambaforge/bin/activate snakemake

    ${params.r} ${projectDir}/bin/proseq_qc.R \\
        --bam ${bam} \\
        --alignment_log ${alignment_log} \\
        --output ${meta.sample}.PROseq_QC_report.txt
    """
}
