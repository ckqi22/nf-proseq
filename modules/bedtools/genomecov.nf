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
    // 链向方向由 params.strandedness 决定（reverse=R1 antisense 标准 PRO-seq / forward=R1 sense）：
    //   reverse: + 基因信号 = 反向比对 read 的 5' 端; single 模式取 -5（R1 5' 端 = 新生 RNA 3' 端 = 活性位点）。
    //   forward: + 基因信号 = 正向比对 read 的 3' 端; single 模式取 -3（R1 3' 端 = 新生 RNA 3' 端 = 活性位点）。
    def strand = params.strandedness?.trim() ?: 'reverse'
    if (strand != 'reverse' && strand != 'forward') { error "params.strandedness must be 'reverse' or 'forward', got: ${strand}" }
    def end_opt  = signal == 'single' ? (strand == 'reverse' ? '-5' : '-3') : ''
    def plus_strand_arg  = (strand == 'reverse') ? '-' : '+'
    def minus_strand_arg = (strand == 'reverse') ? '+' : '-'
    """
    # =========================================================================
    # 链向与端向(gene-strand 约定)
    #   PRO-seq reverse 建库(默认): + 基因信号 = 反向比对 read 的 5'; - 基因信号 = 正向比对 read 的 5'.
    #     forward 建库则对调(由 params.strandedness 控制 -strand 旗标)。
    #   _plus  = + 链基因信号(正值); _minus = - 链基因信号(-scale -1 取负).
    #   signal=single → -5 只取 5' 端(_single 前缀); signal=full → 全长覆盖度(_full 前缀).
    # =========================================================================
    samtools view -H ${bam} | awk '/^@SQ/ {sub(/SN:/, "", \$2); sub(/LN:/, "", \$3); print \$2, \$3}' | LC_COLLATE=C sort -k1,1 > chrom.sizes

    bedtools genomecov -ibam ${bam} ${end_opt} -strand ${plus_strand_arg}  -bg           | LC_COLLATE=C sort -k1,1 -k2,2n > ${meta.sample}_${signal}_plus.bedgraph
    bedtools genomecov -ibam ${bam} ${end_opt} -strand ${minus_strand_arg} -bg -scale -1 | LC_COLLATE=C sort -k1,1 -k2,2n > ${meta.sample}_${signal}_minus.bedgraph

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
