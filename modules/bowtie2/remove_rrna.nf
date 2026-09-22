process REMOVE_RRNA {
    tag "${meta.sample}"

    stageInMode 'symlink'   // 索引目录符号链接 stage（同 align.nf）

    input:
    tuple val(meta), path(reads)       // trimmed reads（SE: [r1]；PE: [r1, r2]）
    tuple val(meta2), path(index_dir)  // rRNA bowtie2 索引父目录（从 *.1.bt2 反推前缀）

    output:
    tuple val(meta), path("*_R{1,2}_rrna_removed.fastq.gz"), emit: clean_reads
    tuple val(meta), path("${meta.sample}_summary_rrna.txt"), emit: rrna_rate

    script:
    // 只按 R1 判定：SE 直接取未比对(--un-gz)；PE 保留「R1 未比对 rRNA」的 read 对（R2 状态忽略）。
    // summary 即 bowtie2 stderr（"overall alignment rate" ≈ rRNA 占比，供 QC）。
    if (meta.single_end) {
        """
        set -euo pipefail
        first=\$(ls ${index_dir}/*.1.bt2 | head -1)
        prefix=\${first%.1.bt2}
        
        bowtie2 \\
            --end-to-end \\
            --very-sensitive \\
            -x \${prefix} \\
            -U ${reads[0]} \\
            --no-unal \\
            --threads 10 \\
            --un-gz ${meta.sample}_R1_rrna_removed.fastq.gz  \\
            > /dev/null 2> ${meta.sample}_summary_rrna.txt
        """
    } else {
        """
        set -euo pipefail
        first=\$(ls ${index_dir}/*.1.bt2 | head -1)
        prefix=\${first%.1.bt2}

        # 只按 R1 判定：保留「R1 未比对 rRNA」的 read 对（flag 4&64）及其 R2 mate（flag 8&128）；
        # R1 比对到 rRNA 的 read 对整体丢弃。samtools view 的 -f/-F 无法表达 OR，用 awk 按 flag 过滤，
        # 再 sort -n 使 mate 相邻 + samtools fastq 回写配对 FASTQ（-0/-s 丢弃不成对残留）。
        bowtie2 \\
            --end-to-end \\
            --very-sensitive \\            
            -x \${prefix} \\
            -1 ${reads[0]} -2 ${reads[1]} \\
            --threads 10 2> ${meta.sample}_summary_rrna.txt \\
            | awk 'BEGIN{OFS="\\t"} /^@/{print; next} {f=\$2+0; if ((int(f/64)%2 && int(f/4)%2) || (int(f/128)%2 && int(f/8)%2)) print}' \\
            | samtools view -bS - \\
            | samtools sort -n - \\
            | samtools fastq -1 ${meta.sample}_R1_rrna_removed.fastq.gz -2 ${meta.sample}_R2_rrna_removed.fastq.gz -0 /dev/null -s /dev/null -N -
        """
    }
}
