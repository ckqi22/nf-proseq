process GTF2SAF {

    container "bio-base:1.0.0"

    input:
    path gtf

    output:
    path 'promoter.saf', emit: promoter_saf
    path 'genebody.saf', emit: genebody_saf
    path 'genebody.bed', emit: genebody_bed
    path 'gene.bed', emit: gene_bed
    path 'plus_genes.bed', emit: plus_genes_bed
    path 'minus_genes.bed', emit: minus_genes_bed

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