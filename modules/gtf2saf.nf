process GTF2SAF {

    container "bio-base:1.0.0"

    input:
    path gtf

    output:
    tuple val("pol2_promoter"), path('promoter.bed'), emit: promoter_bed
    tuple val("pol2_genebody"), path('genebody.bed'), emit: genebody_bed
    tuple val("genebody"), path('genebody_union.saf'), emit: genebody_union_saf
    path 'gene.bed', emit: gene_bed
    // path 'plus_genes.bed', emit: plus_genes_bed     // 未使用，先注释掉
    // path 'minus_genes.bed', emit: minus_genes_bed   // 未使用，先注释掉

    script:
    """
    source /home/ck/miniconda3/bin/activate renv
    Rscript ${projectDir}/bin/gtf2saf.R \\
        --gtf ${gtf} \\
        --tss_upstream ${params.tss.upstream} \\
        --tss_downstream ${params.tss.downstream} \\
        --genebody_offset ${params.tss.genebody_offset} \\
        --outdir ./
    """
}