process SINGLEBASE_COUNT {
    // 单碱基 5' 端定量：对 bedtools genomecov -5 的 bedGraph 按区域求和，
    // 得到每个 gene 在给定区域内的 5' 端（Pol II 活性位点）计数。
    tag "${meta.sample}.${type}"

    input:
    tuple val(meta), path(plus_bg), path(minus_bg)
    tuple val(type), path(region_bed)

    output:
    tuple val(meta), path("${meta.sample}.${type}.counts.txt"), emit: counts

    script:
    """
    source /workplace/hanguojun/mambaforge/bin/activate snakemake

    # 按链拆分区域（第 6 列 = strand），并与 bedGraph 同序排序
    # （sort -k1,1 -k2,2n），保证 bedtools map 两个输入的一致排序。
    awk '\$6=="+"' ${region_bed} | sort -k1,1 -k2,2n > plus.regions.bed
    awk '\$6=="-"' ${region_bed} | sort -k1,1 -k2,2n > minus.regions.bed

    # gene-strand 约定：_plus = + 链基因信号（正值）；_minus = - 链基因信号（负值）。
    # + 基因读 _plus 直接 sum；- 基因读 _minus（负值）sum 后取负回正。
    bedtools map -a plus.regions.bed  -b ${plus_bg}  -c 4 -o sum -null 0 > plus.map
    bedtools map -a minus.regions.bed -b ${minus_bg} -c 4 -o sum -null 0 > minus.map

    # 合并为 gene_id length count（\$4=GeneID，\$3-\$2=区域长度，\$7=map 的 sum 列）。
    awk 'BEGIN{OFS="\\t"} {print \$4, (\$3-\$2),  int(\$7)}' plus.map  > ${meta.sample}.${type}.counts.txt
    awk 'BEGIN{OFS="\\t"} {print \$4, (\$3-\$2), -int(\$7)}' minus.map >> ${meta.sample}.${type}.counts.txt
    sort -k1,1 ${meta.sample}.${type}.counts.txt -o ${meta.sample}.${type}.counts.txt
    """
}
