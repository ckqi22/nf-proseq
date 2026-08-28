process GENEBODY_COUNT {
    tag "${meta.sample}"

    input:
    tuple val(meta), path(plus_bg), path(minus_bg)
    path genebody_bed

    output:
    path "${meta.sample}.genebody.counts.txt", emit: counts

    script:
    """
    source /workplace/hanguojun/mambaforge/bin/activate snakemake

    # 按链拆分 gene body（第 6 列 = strand），并与 bedGraph 同序排序
    # （sort -k1,1 -k2,2n），保证 bedtools map 两个输入的一致排序。
    awk '\$6=="+"' ${genebody_bed} | sort -k1,1 -k2,2n > plus.regions.bed
    awk '\$6=="-"' ${genebody_bed} | sort -k1,1 -k2,2n > minus.regions.bed

    # + 基因：sum plus.bedGraph（read - 链 = + 基因信号，正值）
    bedtools map -a plus.regions.bed -b ${plus_bg} -c 4 -o sum -null 0 > plus.map
    # - 基因：sum minus.bedGraph（read + 链 = - 基因信号，正值）
    bedtools map -a minus.regions.bed -b ${minus_bg} -c 4 -o sum -null 0 > minus.map

    # 合并为 gene_id length count（\$4=GeneID，\$3-\$2=基因体长度，\$7=map 的 sum 列）
    awk 'BEGIN{OFS="\\t"} {print \$4, (\$3-\$2), int(\$7)}' plus.map  > ${meta.sample}.genebody.counts.txt
    awk 'BEGIN{OFS="\\t"} {print \$4, (\$3-\$2), int(\$7)}' minus.map >> ${meta.sample}.genebody.counts.txt
    sort -k1,1 ${meta.sample}.genebody.counts.txt -o ${meta.sample}.genebody.counts.txt
    """
}
