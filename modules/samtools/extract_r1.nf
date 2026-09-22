process EXTRACT_R1 {
    tag "${meta.sample}"

    input:
    tuple val(meta), path(bam)
    val spike_chroms        // spike 染色体名单路径字符串；空串 = 不剔 spike（无 spike-in 时）

    output:
    tuple val(meta), path("${meta.sample}.r1.bam"),           emit: r1_bam
    tuple val(meta), path("${meta.sample}.total_mapped.txt"), emit: total_mapped
    tuple val(meta), path("${meta.sample}.main_mapped.txt"),  emit: main_mapped

    script:
    // total_mapped = 剔 spike 前（含 spike，供 spike-in 占比质控分母）；main_mapped = 剔 spike 后（纯主，供 CPM 分母）。
    // 剔 spike 用 spike_chroms 名单精确匹配（@SQ 头 + body $3），不依赖 spike_ 前缀。
    if (meta.single_end) {
        """
        samtools view -c ${bam} > ${meta.sample}.total_mapped.txt
        if [ -n "${spike_chroms}" ]; then
            samtools view -h ${bam} \\
                | awk 'NR==FNR{sp[\$1]=1; next} /^@SQ/{split(\$2,a,":"); if(a[2] in sp) next} /^@/{print; next} {if(\$3 in sp) next; print}' ${spike_chroms} - \\
                | samtools view -b -o ${meta.sample}.r1.bam -
        else
            cp ${bam} ${meta.sample}.r1.bam
        fi
        samtools view -c ${meta.sample}.r1.bam > ${meta.sample}.main_mapped.txt
        """
    } else {
        """
        samtools view -f 64 -F 4 -c ${bam} > ${meta.sample}.total_mapped.txt
        if [ -n "${spike_chroms}" ]; then
            samtools view -f 64 -F 4 -h ${bam} \\
                | awk 'BEGIN{FS=OFS="\\t"} /^@/{print; next} {\$2=(int(\$2/16)%2)?16:0; print}' \\
                | awk 'NR==FNR{sp[\$1]=1; next} /^@SQ/{split(\$2,a,":"); if(a[2] in sp) next} /^@/{print; next} {if(\$3 in sp) next; print}' ${spike_chroms} - \\
                | samtools view -b -o ${meta.sample}.r1.bam -
        else
            samtools view -f 64 -F 4 -h ${bam} \\
                | awk 'BEGIN{FS=OFS="\\t"} /^@/{print; next} {\$2=(int(\$2/16)%2)?16:0; print}' \\
                | samtools view -b -o ${meta.sample}.r1.bam -
        fi
        samtools view -c ${meta.sample}.r1.bam > ${meta.sample}.main_mapped.txt
        """
    }
}
