process EXTRACT_R1 {
    tag "${meta.sample}"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${meta.sample}.r1.bam"), emit: r1_bam

    script:
    if (meta.single_end) {
        """
        cp ${bam} ${meta.sample}.r1.bam
        """
    } else {
        """
        # PRO-seq: R1 5' 端 = RNA 3' 端 = Pol II 活性位点; R2 是 5' 接头侧、无信号，
        # 任何定量都不应计入 R2。这里抽 read1(-f 64)，再把全部 paired 相关 flag 归零、仅保留 0x10(16, 反链)供下游
        # featureCounts -s 2 判链，输出成"单端"BAM。若不剥 paired flag, featureCounts
        # 单端模式会报 'Paired-end reads were detected in single-end read library'。
        samtools view -f 64 -h ${bam} \\
            | awk 'BEGIN{FS=OFS="\\t"} /^@/{print; next} {\$2=(int(\$2/16)%2)?16:0; print}' \\
            | samtools view -b -o ${meta.sample}.r1.bam -
        """
    }
}
