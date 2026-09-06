process GTF2BED {

    container "bio-base:1.0.0"

    input:
    path gtf
    path rep_gtf

    output:
    path 'tss.bed', emit: tss_bed                                           // 代表 transcript TSS → TSS metagene
    tuple val("pol2_promoter"), path('promoter.bed'), emit: promoter_bed    // 代表 transcript promoter → pol2_count 单碱基 promoter 计数
    tuple val("pol2_genebody"), path('genebody.bed'), emit: genebody_bed    // 代表 transcript genebody → pol2_count 单碱基 genebody 计数
    path 'gene.bed', emit: gene_bed                                         // gene 范围 → SIGNAL_TABLE(信号表 intersect)

    script:
    """
    source /home/ck/miniconda3/bin/activate renv
    Rscript ${projectDir}/bin/gtf2bed.R \\
        --gtf ${gtf} \\
        --rep_gtf ${rep_gtf} \\
        --tss_upstream ${params.tss.upstream} \\
        --tss_downstream ${params.tss.downstream} \\
        --genebody_offset ${params.tss.genebody_offset} \\
        --outdir ./
    """
}
