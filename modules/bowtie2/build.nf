process BOWTIE2_BUILD {
    tag "$fasta"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path('bowtie2'), emit: index

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    mkdir bowtie2
    bowtie2-build \\
    $args \\
    --threads 4 \\
    $fasta \\
    bowtie2/${fasta.baseName}
    """
}
