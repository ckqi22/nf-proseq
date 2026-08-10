process featurecounts_tss {
    tag "${meta.sample}"

    input:
    tuple val(meta), path(plus_bam), path(minus_bam), val(config)

    output:
    path "${meta.sample}.tss_plus.counts.txt",   emit: tss_plus_counts
    path "${meta.sample}.tss_minus.counts.txt",  emit: tss_minus_counts
    path "${meta.sample}.tss_counts.txt",        emit: combined_tss_counts

    script:
    def gtf             = config.gtf
    def upstream        = params.tss?.upstream ?: params.tss_upstream ?: 50
    def downstream      = params.tss?.downstream ?: params.tss_downstream ?: 300
    def feature_counts  = params.feature_counts

    """
    source /workplace/hanguojun/mambaforge/bin/activate snakemake

    # Generate TSS window SAF annotation
    ${params.r} ${projectDir}/bin/generate_tss_regions.R \\
        --gtf ${gtf} \\
        --upstream ${upstream} \\
        --downstream ${downstream} \\
        --output ${meta.sample}.tss_regions.saf

    # PRO-seq strand-specific TSS window quantification (featureCounts -s 2)

    # Plus strand TSS counts
    ${feature_counts} \\
        -s 2 \\
        -F SAF \\
        -g GeneID \\
        -a ${meta.sample}.tss_regions.saf \\
        -o ${meta.sample}.tss_plus.counts.txt \\
        -T 4 \\
        ${plus_bam} 2>&1 | tee -a ${meta.sample}.tss_featureCounts.log

    # Minus strand TSS counts
    ${feature_counts} \\
        -s 2 \\
        -F SAF \\
        -g GeneID \\
        -a ${meta.sample}.tss_regions.saf \\
        -o ${meta.sample}.tss_minus.counts.txt \\
        -T 4 \\
        ${minus_bam} 2>&1 | tee -a ${meta.sample}.tss_featureCounts.log

    # Combine plus and minus strand TSS counts
    ${params.r} ${projectDir}/bin/combine_chain_counts.R \\
        --plus ${meta.sample}.tss_plus.counts.txt \\
        --minus ${meta.sample}.tss_minus.counts.txt \\
        --output ${meta.sample}.tss_counts.txt
    """
}
