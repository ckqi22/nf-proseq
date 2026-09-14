process ALIGN {
    tag "${meta.sample}"

    stageInMode 'symlink'   // 大索引目录用符号链接 stage，避免逐样本拷贝（DB 与 workdir 同盘）

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(index_dir)   // 索引目录（prebuilt 父目录或现建 bowtie2/），从 *.1.bt2 反推前缀


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
    set -euo pipefail

    first=\$(ls ${index_dir}/*.1.bt2 | head -1)
    prefix=\${first%.1.bt2}

    bowtie2 \\
        --end-to-end \\
        --sensitive \\
        ${orientation_args} \\
        -x \${prefix} \\
        ${reads_args} \\
        --threads 10 \\
        2> ${meta.sample}_summary_bowtie2.txt | \\
    samtools view  -@ 2 -bS ${flag_args} -F 4 -q 20 - | \\
    samtools sort  -@ 10 -o ${meta.sample}.bam && \\
    samtools index -@ 1 ${meta.sample}.bam
    """
}