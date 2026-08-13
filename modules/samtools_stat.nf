process SAMTOOLS_STAT {
    tag "${meta.sample}"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(fasta)

    output:
    tuple val(meta), path("*.flagstat"), emit: flagstat
    tuple val(meta), path("*.idxstats"), emit: idxstats
    tuple val(meta), path("*.stats"), emit: stats


    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    samtools \\
        flagstat \\
        --threads 1 \\
        ${bam} \\
        > ${meta.sample}.flagstat
        
    samtools \\
        idxstats \\
        --threads 1 \\
        ${bam} \\
        > ${meta.sample}.idxstats

    samtools \\
        stats \\
        --threads 1 \\
        --reference ${fasta} \\
        ${bam} \\
        > ${meta.sample}.stats
    """    
}