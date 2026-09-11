process SINGLEBASE_COUNT {
    tag "${meta.sample}.${type}"

    input:
    tuple val(meta), path(plus_bigWig), path(minus_bigWig)
    tuple val(type), path(region_bed)

    output:
    tuple val(meta), path("${meta.sample}.${type}.counts.txt"), emit: counts

    script:
    """
    awk '\$6=="+"' ${region_bed} > plus.regions.bed
    awk '\$6=="-"' ${region_bed} > minus.regions.bed

    # 逐区域求和。bigWigAverageOverBed 输出列:
    # 1.name 2.size 3.covered(区域内有覆盖的碱基数, <size) 4.sum(区域内信号和) 5.mean0(sum/size) 6.mean(sum/covered)
    bigWigAverageOverBed ${plus_bigWig}  plus.regions.bed  plus.tab
    bigWigAverageOverBed ${minus_bigWig} minus.regions.bed minus.tab

    # int(x+0.5)=四舍五入，避免 float32 把 15 存成 14.9999 被 int 截成 14
    awk 'BEGIN{OFS="\\t"} {print \$1, \$2, int(\$4+0.5)}'  plus.tab  >  ${meta.sample}.${type}.counts.txt
    awk 'BEGIN{OFS="\\t"} {print \$1, \$2, int(-\$4+0.5)}' minus.tab >> ${meta.sample}.${type}.counts.txt

    sort -k1,1 ${meta.sample}.${type}.counts.txt -o ${meta.sample}.${type}.counts.txt
    """
}


// process SINGLEBASE_COUNT {
//     // 单碱基 5' 端定量：对 bedtools genomecov -5 的 bedGraph 按区域求和，
//     // 得到每个 gene 在给定区域内的 5' 端（Pol II 活性位点）计数。
//     tag "${meta.sample}.${type}"

//     input:
//     tuple val(meta), path(plus_bedgraph), path(minus_bedgraph)
//     tuple val(type), path(region_bed)

//     output:
//     tuple val(meta), path("${meta.sample}.${type}.counts.txt"), emit: counts

//     script:
//     """
//     awk '\$6=="+"' ${region_bed} | sort -k1,1 -k2,2n > plus.regions.bed
//     awk '\$6=="-"' ${region_bed} | sort -k1,1 -k2,2n > minus.regions.bed

//     bedtools intersect -a plus.regions.bed  -b ${plus_bedgraph}  -wa -wb -loj | \
//     awk -F'\\t' 'BEGIN{OFS="\\t"} {g=\$4; L[g]=\$3-\$2; if(\$8 ~ /^[0-9]+\$/){ov=(\$9<\$3?\$9:\$3)-(\$8>\$2?\$8:\$2); s[g]+=\$10*ov}} END{for(g in L) print g, L[g], s[g]+0}' > plus.counts
//     bedtools intersect -a minus.regions.bed -b ${minus_bedgraph} -wa -wb -loj | \
//     awk -F'\\t' 'BEGIN{OFS="\\t"} {g=\$4; L[g]=\$3-\$2; if(\$8 ~ /^[0-9]+\$/){ov=(\$9<\$3?\$9:\$3)-(\$8>\$2?\$8:\$2); s[g]+=-\$10*ov}} END{for(g in L) print g, L[g], s[g]+0}' > minus.counts

//     cat plus.counts minus.counts | sort -k1,1 > ${meta.sample}.${type}.counts.txt

//     """
// }