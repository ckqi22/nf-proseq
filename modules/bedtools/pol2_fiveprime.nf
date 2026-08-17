process POL2_FIVEPRIME {
    tag "${meta.sample}"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${meta.sample}.plus.bedGraph"),   emit: plus_bedgraph
    tuple val(meta), path("${meta.sample}.minus.bedGraph"),  emit: minus_bedgraph
    tuple val(meta), path("${meta.sample}.signed.bedGraph"), emit: signed_bedgraph

    script:
    """
    source /workplace/hanguojun/mambaforge/bin/activate snakemake

    # chrom.sizes from the BAM header (needed by bedtools genomecov -g)
    samtools view -H ${bam} | awk '/^@SQ/ {sub(/SN:/, "", \$2); sub(/LN:/, "", \$3); print \$2, \$3}' > chrom.sizes

    # Collapse each alignment to its 5' end (single base). For a dUTP PRO-seq
    # library this 5' end is the Pol II active site:
    #   + strand -> leftmost (start);   - strand -> rightmost (end)
    bamToBed -i ${bam} \\
        | awk '\$5 > 0 && \$6 == "+" {print \$1, \$2, \$2+1, \$4, \$5, \$6}
               \$5 > 0 && \$6 == "-" {print \$1, \$3-1, \$3, \$4, \$5, \$6}' \\
        | sort -k1,1 -k2,2n > ${meta.sample}.fiveprime.bed

    # Strand-separated coverage (bedGraph); - strand negated (dREG convention)
    bedtools genomecov -bg -i ${meta.sample}.fiveprime.bed -g chrom.sizes -strand + \\
        | sort -k1,1 -k2,2n > ${meta.sample}.plus.bedGraph
    bedtools genomecov -bg -i ${meta.sample}.fiveprime.bed -g chrom.sizes -strand - \\
        | awk '{\$4 = -\$4} 1' \\
        | sort -k1,1 -k2,2n > ${meta.sample}.minus.bedGraph

    # Signed combined bedGraph (dREG input: + positive, - negative)
    cat ${meta.sample}.plus.bedGraph ${meta.sample}.minus.bedGraph \\
        | sort -k1,1 -k2,2n > ${meta.sample}.signed.bedGraph
    """
}
