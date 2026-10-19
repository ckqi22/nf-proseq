process ALIGN {
    tag "${meta.sample}"

    stageInMode 'symlink'   // 大索引目录用符号链接 stage，避免逐样本拷贝（DB 与 workdir 同盘）

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(index_dir)   // 索引目录（prebuilt 父目录或现建 bowtie2/），从 *.1.bt2 反推前缀


    output:
    tuple val(meta), path("${meta.sample}.bam"), emit: bam
    tuple val(meta), path("${meta.sample}.bam.bai"), emit: bai
    tuple val(meta), path("${meta.sample}_unmapped.bam"), emit: unmapped_bam
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
    // 未比对 reads 单独存：SE=整条未比对(-f 4)；PE=只取 R1 未比对(-f 68=64|4，不要求 R2 状态)。
    def unmapped_flag = meta.single_end ? "4" : "68"
    """
    set -euo pipefail

    first=\$(ls ${index_dir}/*.1.bt2 | head -1)
    prefix=\${first%.1.bt2}

    # bowtie2 输出 SAM(含未比对)到中间文件，再从同一份 SAM 拆两条：
    #   主 bam = mapped + MAPQ>=20(过滤条件与改前一致，主 bam 字节不变);
    #   _unmapped.bam = R1 未比对 reads(SE 即整条未比对)，供排查/重比对。
    bowtie2 \\
        --end-to-end \\
        --very-sensitive \\
        ${orientation_args} \\
        -x \${prefix} \\
        ${reads_args} \\
        --threads ${task.cpus} \\
        -S ${meta.sample}.sam \\
        2> ${meta.sample}_summary_bowtie2.txt

    # 比对率 QC：overall alignment rate 低于阈值直接报错终止（失败快速，避免后续 sort 白跑）
    awk -v t="${params.threshold.mapping_rate}" '
        /overall alignment rate/ {
            gsub("%", "", \$1); rate = \$1 + 0;
            print "  [QC] ${meta.sample}: overall alignment rate = " rate "%";
            if (rate < t * 100) {
                print "ERROR: ${meta.sample}: alignment rate " rate "% below threshold " (t * 100) "%" > "/dev/stderr";
                exit 1;
            }
            found = 1;
        }
        END { if (!found) { print "ERROR: overall alignment rate not found in bowtie2 summary" > "/dev/stderr"; exit 1 } }
    ' ${meta.sample}_summary_bowtie2.txt

    samtools view -@ 2 -bS ${flag_args} -F 4 -q 20 ${meta.sample}.sam \\
        | samtools sort -@ ${task.cpus} -o ${meta.sample}.bam
    samtools view -@ 2 -bS -f ${unmapped_flag} ${meta.sample}.sam \\
        | samtools sort -@ ${task.cpus} -o ${meta.sample}_unmapped.bam
    samtools index -@ 1 ${meta.sample}.bam
    rm -f ${meta.sample}.sam
    """
}