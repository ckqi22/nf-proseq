process featurecounts_gb {
    tag "${meta.sample}"

    input:
    tuple val(meta), path(plus_bam), path(minus_bam), val(config)

    output:
    path "${meta.sample}.gene_body_plus.counts.txt",   emit: plus_counts
    path "${meta.sample}.gene_body_minus.counts.txt",  emit: minus_counts
    path "${meta.sample}.gene_body_counts.txt",        emit: combined_counts

    script:
    def gtf            = config.gtf
    def offset         = params.tss?.gene_body_offset ?: params.gene_body_offset ?: 301
    def feature_counts = params.feature_counts

    """
    source /workplace/hanguojun/mambaforge/bin/activate snakemake

    # PRO-seq strand-specific gene body quantification (featureCounts -s 2 reverse-stranded)
    # Gene body: TSS+${offset} to TES (full gene minus the TSS-proximal pause window)

    # Plus strand BAM
    ${feature_counts} \\
        -s 2 \\
        -t exon \\
        -g gene_id \\
        -a ${gtf} \\
        -o ${meta.sample}.gene_body_plus.counts.txt \\
        -T 4 \\
        ${plus_bam} 2>&1 | tee -a ${meta.sample}.gene_body_featureCounts.log

    # Minus strand BAM
    ${feature_counts} \\
        -s 2 \\
        -t exon \\
        -g gene_id \\
        -a ${gtf} \\
        -o ${meta.sample}.gene_body_minus.counts.txt \\
        -T 4 \\
        ${minus_bam} 2>&1 | tee -a ${meta.sample}.gene_body_featureCounts.log

    # Combine plus and minus strand counts
    ${params.r} ${projectDir}/bin/combine_chain_counts.R \\
        --plus ${meta.sample}.gene_body_plus.counts.txt \\
        --minus ${meta.sample}.gene_body_minus.counts.txt \\
        --output ${meta.sample}.gene_body_counts.txt
    """
}
