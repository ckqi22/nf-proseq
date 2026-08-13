process FEATURECOUNTS {
    tag "${meta.sample}"

    input:
    tuple val(meta), path(bams)
    path annotation
    val type

    output:
    path "${meta.sample}.proseq_counts.txt", emit: counts
    path "${meta.sample}.proseq_counts.txt.summary", emit: summary

    script:
    def paired_end = meta.single_end ? '' : '-p'

    def strandedness = 0
    if (params.strandedness == 'forward') {
        strandedness = 1
    } else if (params.strandedness == 'reverse') {
        strandedness = 2
    }

    def gtf             = config.gtf
    def feature_counts  = params.feature_counts
    def upstream        = params.tss?.upstream ?: params.tss_upstream ?: 50
    def downstream      = params.tss?.downstream ?: params.tss_downstream ?: 300
    def offset          = params.tss?.gene_body_offset ?: params.gene_body_offset ?: 301


    """
    source /workplace/hanguojun/mambaforge/bin/activate snakemake

    ${params.feature_counts} \\
        ${paired_end} \\
        -T 2 \\
        -a ${annotation} \\
        -t ${type} \\
        -g gene_id \\
        -s ${strandedness} \\
        -o ${prefix}.featureCounts.txt \\
        ${bams.join(' ')}
















    # === Step 1: Generate both SAF files (one GTF parse) ===

    ${params.r} ${projectDir}/bin/generate_proseq_saf.R \\
        --gtf ${gtf} \\
        --upstream ${upstream} \\
        --downstream ${downstream} \\
        --offset ${offset} \\
        --tss_out ${meta.sample}.tss.saf \\
        --gb_out ${meta.sample}.gene_body.saf

    # === Step 2: Run featureCounts on three regions ===

    # ① Full gene (GTF, -t gene)
    ${feature_counts} \\
        -s 2 -p -t gene -g gene_id \\
        -a ${gtf} \\
        -o ${meta.sample}.full_gene.txt \\
        -T 4 \\
        ${bam} 2>&1 | tee ${meta.sample}.full_gene.log

    # ② TSS window (SAF)
    ${feature_counts} \\
        -s 2 -p -F SAF -g GeneID \\
        -a ${meta.sample}.tss.saf \\
        -o ${meta.sample}.tss.txt \\
        -T 4 \\
        ${bam} 2>&1 | tee ${meta.sample}.tss.log

    # ③ Gene body (SAF)
    ${feature_counts} \\
        -s 2 -p -F SAF -g GeneID \\
        -a ${meta.sample}.gene_body.saf \\
        -o ${meta.sample}.gene_body.txt \\
        -T 4 \\
        ${bam} 2>&1 | tee ${meta.sample}.gene_body.log

    # === Step 3: Merge into one file ===
    ${params.r} -e "
        full  <- read.delim('${meta.sample}.full_gene.txt',  header=TRUE, stringsAsFactors=FALSE, comment.char='#')
        tss   <- read.delim('${meta.sample}.tss.txt',        header=TRUE, stringsAsFactors=FALSE, comment.char='#')
        gb    <- read.delim('${meta.sample}.gene_body.txt',  header=TRUE, stringsAsFactors=FALSE, comment.char='#')

        # Align by gene_id
        col_full <- tail(names(full), 1)
        col_tss  <- tail(names(tss),  1)
        col_gb   <- tail(names(gb),   1)

        result <- full[, c('Geneid', 'Chr', 'Start', 'End', 'Strand', 'Length'), drop=FALSE]
        names(result)[1] <- 'gene_id'

        result[[paste0('${meta.sample}', '_full_gene')]]  <- full[[col_full]] [match(result\$gene_id, full\$Geneid)]
        result[[paste0('${meta.sample}', '_tss')]]        <- tss[[col_tss]]  [match(result\$gene_id, tss\$Geneid)]
        result[[paste0('${meta.sample}', '_gene_body')]]  <- gb[[col_gb]]   [match(result\$gene_id, gb\$Geneid)]

        result[is.na(result)] <- 0

        write.table(result, file='${meta.sample}.proseq_counts.txt',
                    sep='\\t', quote=FALSE, row.names=FALSE)
        message('[merge] ', nrow(result), ' genes with 3 region counts written')
    "
    """
}
