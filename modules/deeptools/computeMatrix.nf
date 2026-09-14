process COMPUTEMATRIX {
    tag "${meta.sample}_${label}"

    conda '/workplace/hanguojun/mambaforge/envs/deeptools'
    cpus 4

    input:
    tuple val(meta), path(plus_bw), path(minus_bw)
    path  tss_bed
    val   upstream        // TSS 上游窗口（bp）
    val   downstream      // TSS 下游窗口（bp）
    val   bin_size        // bin 大小（bp）
    val   label           // 'profile' | 'heatmap'：区分两条矩阵链，避免输出名相撞

    output:
    tuple val(meta), path("${meta.sample}_${label}_TSS_merged_matrix.gz"), emit: matrix

    script:
    // 正/负链基因分开算（各用对应链的 bigWig），再 rbind 合并为单 group 矩阵。
    // 兼容单样本（单个 Path）与多样本叠加（List<Path>）：统一归一化为列表；
    // 样本标签取 meta.sample_names（叠加时由 metagene.nf / metagene_group.nf 写入），
    // 单样本退化为 [meta.sample]。
    def plus_list  = (plus_bw instanceof List) ? plus_bw : [plus_bw]
    def minus_list = (minus_bw instanceof List) ? minus_bw : [minus_bw]
    def names      = (meta.sample_names instanceof List) ? meta.sample_names : [meta.sample]
    """
    awk '\$6=="+"' ${tss_bed} > regions_plus.bed
    awk '\$6=="-"' ${tss_bed} > regions_minus.bed

    computeMatrix reference-point \\
        -R regions_plus.bed \\
        -S ${plus_list.join(' ')} \\
        --referencePoint TSS \\
        --upstream ${upstream} \\
        --downstream ${downstream} \\
        --binSize ${bin_size} \\
        --missingDataAsZero \\
        --samplesLabel ${names.join(' ')} \\
        -p ${task.cpus} \\
        -o ${meta.sample}_plus_${label}_TSS_matrix.gz

    computeMatrix reference-point \\
        -R regions_minus.bed \\
        -S ${minus_list.join(' ')} \\
        --referencePoint TSS \\
        --upstream ${upstream} \\
        --downstream ${downstream} \\
        --binSize ${bin_size} \\
        --scale -1 \\
        --missingDataAsZero \\
        --samplesLabel ${names.join(' ')} \\
        -p ${task.cpus} \\
        -o ${meta.sample}_minus_${label}_TSS_matrix.gz

    computeMatrixOperations rbind \\
        -m ${meta.sample}_plus_${label}_TSS_matrix.gz ${meta.sample}_minus_${label}_TSS_matrix.gz \\
        -o ${meta.sample}_${label}_TSS_merged_matrix.gz
    """
}

// 说明：
// - tss.bed（bin/gtf2bed.R 产出）为 BED6 无表头：Start=tss-1, End=tss，第 6 列 strand。
//   deepTools reference-point TSS：+ 链以 region start、- 链以 region end 为参考点，
//   且 - 链矩阵自动反转为 5'→3'（heatmapper.py 中 cov[::-1]），故两链矩阵方向一致
//   （TSS 在左、沿基因 5'→3'）。上游/下游自动延伸出 region 边界。
// - minus_cpm bigWig 由 genomecov `-scale -1` 写成负值，minus 矩阵 --scale -1 翻正；
//   NaN×(-1) 仍为 NaN，被 --missingDataAsZero 填 0，与 R 版 abs-逐位置-再 bin 数学等价。
// - 矩阵列数 = (upstream + downstream) / binSize，如 ±1000 / 10 → 200 列。
// - rbind 需 deepTools >= 3.1.2（v3.1.1 有 group_labels 未更新的 bug）。
// - 多样本叠加：metagene.nf / metagene_group.nf 把「全部样本/组的 bigWig」收集成 List
//   后传 -S，两侧 --samplesLabel 同序，rbind 按 region 行堆叠、保留 N 个 sample 列 → plotProfile 一图多线。
