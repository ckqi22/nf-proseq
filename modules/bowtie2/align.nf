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
    // PRO-seq 链向：R1 落在反义链(reverse)、R2 落在正义链(forward)，即 mate1 反链、mate2 正链 = --rf
    def orientation_args = meta.single_end ? "" : "--rf"
    """
    source /workplace/hanguojun/mambaforge/bin/activate snakemake

    bowtie2 \\
        --end-to-end \\
        --sensitive \\
        ${orientation_args} \\
        -x ${index} \\
        ${reads_args} \\
        --threads 10 \\
        2> >(tee ${meta.sample}_summary_bowtie2.txt) | \\
    samtools view  -@ 2 -bS ${flag_args} -q 20 - | \\
    samtools sort  -@ 10 -o ${meta.sample}.bam && \\
    samtools index -@ 1 ${meta.sample}.bam
    """
}