process ALIGN {
    tag "${meta.sample}"

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), val(index)


    output:
    tuple val(meta), path("${meta.sample}.bam"), emit: bam
    tuple val(meta), path("${meta.sample}.bam.bai"), emit: bai
    tuple val(meta), path("${meta.sample}_summary_bowtie2.txt"), emit: alignRate


    script:
    def reads_args = ""
    def flag_args = ""
    if (meta.single_end) {
        reads_args = "-U ${reads}"
    } else {
        reads_args = "-1 ${reads[0]} -2 ${reads[1]}"
        flag_args = "-f 2"
    }
    def orientation_args = meta.single_end ? "" : "--fr"
    """
    bowtie2 \\
        --end-to-end \\
        --sensitive \\
        ${orientation_args} \\
        -x ${index} \\
        ${reads_args} \\
        --threads 10 \\
        2> >(tee ${meta.sample}_summary_bowtie2.txt) | \\
    samtools view  -@ 2 -bS ${flag_args} -F 4 -q 20 - | \\
    samtools sort  -@ 10 -o ${meta.sample}.bam && \\
    samtools index -@ 1 ${meta.sample}.bam
    """
}