process GENOMECOV {
    tag "${meta.sample}_${signal}"

    input:
    tuple val(meta), path(bam)
    val signal   // 'single'（5' 端单碱基）| 'full'（全长覆盖度）——文件名前缀直接取该值

    output:
    tuple val(meta), path("${meta.sample}_${signal}_plus.bedgraph"),      path("${meta.sample}_${signal}_minus.bedgraph"),      emit: bedgraph
    tuple val(meta), path("${meta.sample}_${signal}_plus.bigWig"),        path("${meta.sample}_${signal}_minus.bigWig"),        emit: bigwig
    tuple val(meta), path("${meta.sample}_${signal}_plus_cpm.bedgraph"),  path("${meta.sample}_${signal}_minus_cpm.bedgraph"),  emit: bedgraph_cpm
    tuple val(meta), path("${meta.sample}_${signal}_plus_cpm.bigWig"),    path("${meta.sample}_${signal}_minus_cpm.bigWig"),    emit: bigwig_cpm
    tuple val(meta), path("${meta.sample}_${signal}_plus_spike.bigWig"),  path("${meta.sample}_${signal}_minus_spike.bigWig"),  emit: bigwig_spike, optional: true

    script:
    // 由 signal 派生：single → -5（5' 端）；full → 全长覆盖度。文件名前缀 _single/_full 直接取 signal。
    // 进程只算一种信号，门控在调用侧（pol2_count 按 mode 决定调 GENOMECOV 还是 GENOMECOV_FULL 别名）。
    def five = signal == 'single' ? '-5' : ''
    """
    # =========================================================================
    # 链向与端向(gene-strand 约定)
    #   PRO-seq reverse 建库：+ 基因信号 = 反向比对 read 的 5'; - 基因信号 = 正向比对 read 的 5'.
    #   _plus  = + 链基因信号(-strand -, 正值); _minus = - 链基因信号(-strand +, -scale -1 取负).
    #   signal=single → -5 只取 5' 端(_single 前缀); signal=full → 全长覆盖度(_full 前缀).
    # spike-in: combined BAM 含 spike_* 染色体; 下游只关心主基因组, chrom.sizes 与各
    #   bedGraph 一律剔除 spike_* 染色体，使 bigWig/bedGraph 为纯主基因组。
    # =========================================================================
    samtools view -H ${bam} | awk '/^@SQ/ {sub(/SN:/, "", \$2); sub(/LN:/, "", \$3); print \$2, \$3}' | awk '\$1 !~ /^spike_/' | sort -k1,1 > chrom.sizes

    bedtools genomecov -ibam ${bam} ${five} -strand - -bg           | awk '\$1 !~ /^spike_/' | sort -k1,1 -k2,2n > ${meta.sample}_${signal}_plus.bedgraph
    bedtools genomecov -ibam ${bam} ${five} -strand + -bg -scale -1 | awk '\$1 !~ /^spike_/' | sort -k1,1 -k2,2n > ${meta.sample}_${signal}_minus.bedgraph

    awk -v s="${meta.scale_cpm}" 'BEGIN{OFS="\\t"}{\$4=\$4*s; print}' ${meta.sample}_${signal}_plus.bedgraph  > ${meta.sample}_${signal}_plus_cpm.bedgraph
    awk -v s="${meta.scale_cpm}" 'BEGIN{OFS="\\t"}{\$4=\$4*s; print}' ${meta.sample}_${signal}_minus.bedgraph > ${meta.sample}_${signal}_minus_cpm.bedgraph

    bedGraphToBigWig ${meta.sample}_${signal}_plus.bedgraph       chrom.sizes ${meta.sample}_${signal}_plus.bigWig
    bedGraphToBigWig ${meta.sample}_${signal}_minus.bedgraph      chrom.sizes ${meta.sample}_${signal}_minus.bigWig
    bedGraphToBigWig ${meta.sample}_${signal}_plus_cpm.bedgraph   chrom.sizes ${meta.sample}_${signal}_plus_cpm.bigWig
    bedGraphToBigWig ${meta.sample}_${signal}_minus_cpm.bedgraph  chrom.sizes ${meta.sample}_${signal}_minus_cpm.bigWig

    # spike 版：同一 raw bedgraph 换 scale(1e6/spike_count); 开 spike 才有(scale_spike != null)
    ${meta.scale_spike != null ? """
    awk -v s="${meta.scale_spike}" 'BEGIN{OFS="\\t"}{\$4=\$4*s; print}' ${meta.sample}_${signal}_plus.bedgraph  > ${meta.sample}_${signal}_plus_spike.bedgraph
    awk -v s="${meta.scale_spike}" 'BEGIN{OFS="\\t"}{\$4=\$4*s; print}' ${meta.sample}_${signal}_minus.bedgraph > ${meta.sample}_${signal}_minus_spike.bedgraph
    bedGraphToBigWig ${meta.sample}_${signal}_plus_spike.bedgraph  chrom.sizes ${meta.sample}_${signal}_plus_spike.bigWig
    bedGraphToBigWig ${meta.sample}_${signal}_minus_spike.bedgraph chrom.sizes ${meta.sample}_${signal}_minus_spike.bigWig
    """ : ''}
    """
}


// bamCoverage \
//  -b HS0_rep2.bam \
//  -o HS0_rep2_plus_deeptools.bedgraph \
//  --binSize 1 \
//  --filterRNAstrand forward \
//  --Offset 1 \
//  --outFileFormat bedgraph \
// 	--numberOfProcessors 10 \
//  --skipNAs
// 与 bedtools genomecov -ibam ${bam} -5 -strand - -bg 一致（bamCoverage --filterRNAstrand forward = FLAG16 反向比对 = -strand -，实测正确）
