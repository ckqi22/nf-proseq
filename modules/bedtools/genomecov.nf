process GENOMECOV {
    tag "${meta.sample}"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${meta.sample}_plus.bedgraph"), path("${meta.sample}_minus.bedgraph"), emit: bedgraph
    tuple val(meta), path("${meta.sample}_plus.bigWig"), path("${meta.sample}_minus.bigWig"), emit: bigwig
    tuple val(meta), path("${meta.sample}_plus_cpm.bedgraph"), path("${meta.sample}_minus_cpm.bedgraph"), emit: bedgraph_cpm
    tuple val(meta), path("${meta.sample}_plus_cpm.bigWig"), path("${meta.sample}_minus_cpm.bigWig"), emit: bigwig_cpm
    tuple val(meta), path("${meta.sample}_forward_cpm.bigWig"), path("${meta.sample}_reverse_cpm.bigWig"), emit: bigwig_coverage_cpm

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
    // PE 时只保留 read1（flag 0x40=64，R2 是 5' 接头侧无信号）；SE 直接使用 BAM。
    // The input bedGraph file must be sorted
    def sig_bam = meta.single_end ? "${bam}" : "r1.bam"
    def extract = meta.single_end ? "" : "samtools view -f 64 -F 4 -b ${bam} -o r1.bam"

    """
    samtools view -H ${bam} | awk '/^@SQ/ {sub(/SN:/, "", \$2); sub(/LN:/, "", \$3); print \$2, \$3}' | sort -k1,1 > chrom.sizes

    ${extract}

    total=\$(samtools view -c ${sig_bam})
    scale=\$(awk -v t="\$total" 'BEGIN{printf "%.12g", 1000000.0/t}')

    bedtools genomecov -ibam ${sig_bam} -5 -strand - -bg           | sort -k1,1 -k2,2n > ${meta.sample}_plus.bedgraph
    bedtools genomecov -ibam ${sig_bam} -5 -strand + -bg -scale -1 | sort -k1,1 -k2,2n > ${meta.sample}_minus.bedgraph

    awk -v s="\$scale" 'BEGIN{OFS="\\t"}{\$4=\$4*s; print}' ${meta.sample}_plus.bedgraph  > ${meta.sample}_plus_cpm.bedgraph
    awk -v s="\$scale" 'BEGIN{OFS="\\t"}{\$4=\$4*s; print}' ${meta.sample}_minus.bedgraph > ${meta.sample}_minus_cpm.bedgraph

    bedGraphToBigWig ${meta.sample}_plus.bedgraph       chrom.sizes ${meta.sample}_plus.bigWig
    bedGraphToBigWig ${meta.sample}_minus.bedgraph      chrom.sizes ${meta.sample}_minus.bigWig
    bedGraphToBigWig ${meta.sample}_plus_cpm.bedgraph   chrom.sizes ${meta.sample}_plus_cpm.bigWig
    bedGraphToBigWig ${meta.sample}_minus_cpm.bedgraph  chrom.sizes ${meta.sample}_minus_cpm.bigWig

    # read全长覆盖度
    # forward = 正链基因信号(reverse read, 正值); reverse = 负链基因信号(forward read, 负值)
    bedtools genomecov -ibam ${sig_bam} -strand - -bg           | sort -k1,1 -k2,2n > ${meta.sample}_forward.bedgraph
    bedtools genomecov -ibam ${sig_bam} -strand + -bg -scale -1 | sort -k1,1 -k2,2n > ${meta.sample}_reverse.bedgraph

    awk -v s="\$scale" 'BEGIN{OFS="\t"}{\$4=\$4*s; print}' ${meta.sample}_forward.bedgraph > ${meta.sample}_forward_cpm.bedgraph
    awk -v s="\$scale" 'BEGIN{OFS="\t"}{\$4=\$4*s; print}' ${meta.sample}_reverse.bedgraph > ${meta.sample}_reverse_cpm.bedgraph

    bedGraphToBigWig ${meta.sample}_forward_cpm.bedgraph chrom.sizes ${meta.sample}_forward_cpm.bigWig
    bedGraphToBigWig ${meta.sample}_reverse_cpm.bedgraph chrom.sizes ${meta.sample}_reverse_cpm.bigWig
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
// 与bedtools genomecov -ibam ${sig_bam} -5 -strand - -bg一致