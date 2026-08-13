process FASTQC {
    tag "${meta.sample}"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*fastqc.zip"), emit: zip
    tuple val(meta), path("*fastqc.html"), emit: html

    script:
    def args = task.ext.args ?: ""
    """
    fastqc \\
        ${args} \\
        --threads 4 \\
        ${reads}
    """
}