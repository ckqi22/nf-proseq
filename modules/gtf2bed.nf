process GTF2BED {

    container "bio-base:1.0.0"

    input:
    path gtf

    output:
    path 'tss.bed', emit: tss_bed                                           // 各 protein_coding transcript TSS → TSS metagene
    tuple val("pol2_promoter"), path('promoter.bed'), emit: promoter_bed    // 各 transcript promoter → pol2_count 单碱基 promoter 计数
    tuple val("pol2_genebody"), path('genebody.bed'), emit: genebody_bed    // 各 transcript genebody → pol2_count 单碱基 genebody 计数
    path 'tx2gene.tsv', emit: tx2gene                                       // transcript_id → gene_id 映射 → merge 加 gene_id 列

    script:
    """
    source /home/ck/miniconda3/bin/activate renv
    Rscript ${projectDir}/bin/gtf2bed.R \\
        --gtf ${gtf} \\
        --tss_upstream ${params.tss.upstream} \\
        --tss_downstream ${params.tss.downstream} \\
        --genebody_offset ${params.tss.genebody_offset} \\
        --outdir ./
    """
}
