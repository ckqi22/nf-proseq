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
    // 正/负链基因分开算（各用对应链的 bigWig），再 rbind 合并为单 group 矩阵：
    """
    awk '\$6=="+"' ${tss_bed} > regions_plus.bed
    awk '\$6=="-"' ${tss_bed} > regions_minus.bed

    computeMatrix reference-point \\
        -R regions_plus.bed \\
        -S ${plus_bw} \\
        --referencePoint TSS \\
        --upstream ${upstream} \\
        --downstream ${downstream} \\
        --binSize ${bin_size} \\
        --missingDataAsZero \\
        --samplesLabel ${meta.sample} \\
        -p ${task.cpus} \\
        -o ${meta.sample}_plus_${label}_TSS_matrix.gz

    computeMatrix reference-point \\
        -R regions_minus.bed \\
        -S ${minus_bw} \\
        --referencePoint TSS \\
        --upstream ${upstream} \\
        --downstream ${downstream} \\
        --binSize ${bin_size} \\
        --scale -1 \\
        --missingDataAsZero \\
        --samplesLabel ${meta.sample} \\
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
