process GENOMECOV {
    tag "${meta.sample}"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${meta.sample}_plus.bedgraph"), path("${meta.sample}_minus.bedgraph"), emit: bedgraph
    tuple val(meta), path("${meta.sample}_plus.bigWig"), path("${meta.sample}_minus.bigWig"), emit: bigwig

    script:
    // =========================================================================
    // 链向与端向（gene-strand 约定）
    // -------------------------------------------------------------------------
    // PRO-seq 建库 reverse：+ 基因信号落在反向比对 read 的 5'、- 基因信号落在正向比对 read 的 5'。
    // gene-strand 命名（按基因链）：
    //   _plus  = + 链基因信号（reverse，-strand -） => 正值
    //   _minus = - 链基因信号（forward，-strand +） => 取负（-scale -1）
    // 两轨一正一负，加载进 IGV 即 gene-strand signed 视图（同官方 TrackTx/Mahat）。
    // 计数时按基因链取对应轨（pol2_count），负值轨取负回正。
    // PE 时只保留 read1（flag 0x40=64，R2 是 5' 接头侧无信号）；SE 直接使用 BAM。
    def sig_bam = meta.single_end ? "${bam}" : "r1.bam"
    def extract = meta.single_end ? "" : "samtools view -f 64 -F 4 -b ${bam} -o r1.bam"

    """
    source /workplace/hanguojun/mambaforge/bin/activate snakemake

    samtools view -H ${bam} | awk '/^@SQ/ {sub(/SN:/, "", \$2); sub(/LN:/, "", \$3); print \$2, \$3}' | sort -k1,1 > chrom.sizes

    ${extract}

    bedtools genomecov -ibam ${sig_bam} -5 -strand - -bg           | sort -k1,1 -k2,2n > ${meta.sample}_plus.bedgraph
    bedtools genomecov -ibam ${sig_bam} -5 -strand + -bg -scale -1 | sort -k1,1 -k2,2n > ${meta.sample}_minus.bedgraph

    bedGraphToBigWig ${meta.sample}_plus.bedgraph  chrom.sizes ${meta.sample}_plus.bigWig
    bedGraphToBigWig ${meta.sample}_minus.bedgraph chrom.sizes ${meta.sample}_minus.bigWig
    """
}
