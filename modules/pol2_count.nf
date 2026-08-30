process POL2_COUNT {
    tag "${meta.sample}"

    input:
    tuple val(meta), path(plus_bedGraph), path(minus_bedGraph)
    path gene_bed

    output:
    path "${meta.sample}.pol2.counts.txt", emit: counts

    script:
    """
    awk '\$6=="+"' ${gene_bed} | sort -k1,1 -k2,2n > plus.regions.bed
    awk '\$6=="-"' ${gene_bed} | sort -k1,1 -k2,2n > minus.regions.bed

    # gene-strand: _plus = + 链基因信号(正); _minus = - 链基因信号(取负)
    # + 链基因信号 = plus_bedGraph(正值) → 直接 sum。
    bedtools map \\
    -a plus.regions.bed \\
    -b ${plus_bedGraph} \\
    -c 4 \\
    -o sum \\
    -null 0 > plus.map

    # - 链基因信号 = minus_bedGraph(负值) → sum 后取负回正。
    bedtools map \\
    -a minus.regions.bed \\
    -b ${minus_bedGraph} \\
    -c 4 \\
    -o sum \\
    -null 0 > minus.map

    awk 'BEGIN{OFS="\\t"} {print \$4, (\$3-\$2),  int(\$7)}' plus.map  > ${meta.sample}.pol2.counts.txt
    awk 'BEGIN{OFS="\\t"} {print \$4, (\$3-\$2), -int(\$7)}' minus.map >> ${meta.sample}.pol2.counts.txt
    
    awk 'NR==FNR{order[\$4]=NR; next} {print order[$1]"\t"$0}' \
    ${gene_bed} ${meta.sample}.pol2.counts.txt \
    | sort -k1,1n \
    | cut -f2- \
    > ${meta.sample}.pol2.counts.sorted.txt
    """
}
