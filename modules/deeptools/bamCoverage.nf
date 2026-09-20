process BAMCOVERAGE {
    tag "${meta.sample}"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${meta.sample}_plus.bedgraph"), path("${meta.sample}_minus.bedgraph"), emit: bedgraph
    tuple val(meta), path("${meta.sample}_plus.bigWig"), path("${meta.sample}_minus.bigWig"), emit: bigwig
    tuple val(meta), path("${meta.sample}_plus_cpm.bedgraph"), path("${meta.sample}_minus_cpm.bedgraph"), emit: bedgraph_cpm
    tuple val(meta), path("${meta.sample}_plus_cpm.bigWig"), path("${meta.sample}_minus_cpm.bigWig"), emit: bigwig_cpm
    tuple val(meta), path("${meta.sample}_full_plus_cpm.bigWig"), path("${meta.sample}_full_minus_cpm.bigWig"), emit: bigwig_coverage_cpm
    tuple val(meta), path("${meta.sample}_full_plus.bigWig"), path("${meta.sample}_full_minus.bigWig"), emit: bigwig_full

    script:
    // =========================================================================
    // 链向与端向（gene-strand 约定）
    // -------------------------------------------------------------------------
    // PRO-seq 建库 reverse：+ 基因信号落在反向比对 read 的 5'、- 基因信号落在正向比对 read 的 5'。
    // gene-strand 命名（按基因链）：
    //   _plus  = + 链基因信号（reverse，-strand -） => 正值
    //   _minus = - 链基因信号（forward，-strand +） => 取负（-scale -1）
    // 两轨一正一负，加载进 IGV 即 gene-strand signed 视图（同官方 TrackTx/Mahat）。
    // 计数时按基因链取对应轨（singlebase_count），负值轨取负回正。
    // 输入 bam 已是 read1 单端化 BAM（align_bowtie2 的 EXTRACT_R1 产出；SE=原 bam 拷贝），
    // 此处不再自行抽 read1。
    // The input bedGraph file must be sorted
    // CPM scale 由 pol2_count 依据 denom（spike_count || total_mapped）算好，经 meta.scale_cpm 传入
    // （此模块未接线；启用时用 --scaleFactor ${meta.scale_cpm}）。

    """
    # spike-in：combined BAM 含 spike_* 染色体（spikein_concat.nf 统一加 spike_ 前缀）。
    # 下游只关心主基因组，chrom.sizes 与各 bedGraph 一律剔除 spike_* 染色体，使 bigWig/bedGraph 为纯主基因组。
    samtools view -H ${bam} | awk '/^@SQ/ {sub(/SN:/, "", \$2); sub(/LN:/, "", \$3); print \$2, \$3}' | awk '\$1 !~ /^spike_/' | sort -k1,1 > chrom.sizes

    # single base coverage
    bamCoverage \\
        --bam ${bam} \\
        --outFileName ${meta.sample}_plus.bedgraph \\
        --outFileFormat bedgraph \\
        --scaleFactor 1 \\
        --Offset 1 \\
        --filterRNAstrand forward \\
        --binSize 1 \\
        --numberOfProcessors 10 \\
        --normalizeUsing None \\
        --skipNAs

    bamCoverage \\
        --bam ${bam} \\
        --outFileName ${meta.sample}_minus.bedgraph \\
        --outFileFormat bedgraph \\
        --scaleFactor -1 \\
        --Offset 1 \\
        --filterRNAstrand reverse \\
        --binSize 1 \\
        --numberOfProcessors 10 \\
        --normalizeUsing None \\
        --skipNAs

    # full read coverage
    
    """
}
