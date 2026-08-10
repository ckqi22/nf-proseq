process metagene_tss {
    tag "metagene_tss"

    input:
    path plus_bam_list
    path minus_bam_list
    val config

    output:
    path "TSS_Metagene_Plus.pdf",   emit: tss_plus_pdf
    path "TSS_Metagene_Minus.pdf",  emit: tss_minus_pdf

    script:
    def window = params.tss_metagene_window ?: 3000

    """
    source /workplace/hanguojun/mambaforge/bin/activate snakemake

    # TSS metagene plot (strand-specific) using deepTools
    # Plus strand BAMs on plus-strand genes; Minus strand BAMs on minus-strand genes

    # Convert BAM lists to space-separated strings
    PLUS_BAMS=\$(cat ${plus_bam_list} | tr '\\n' ' ')
    MINUS_BAMS=\$(cat ${minus_bam_list} | tr '\\n' ' ')

    # Plus strand: computeMatrix and plot
    computeMatrix scale-regions \\
        -S \$PLUS_BAMS \\
        -R ${config.gtf} \\
        --beforeRegionStartLength ${window} \\
        --regionBodyLength 1000 \\
        --afterRegionStartLength ${window} \\
        --skipZeros \\
        -o TSS_Plus_matrix.gz \\
        -p 8

    plotProfile \\
        -m TSS_Plus_matrix.gz \\
        -o TSS_Metagene_Plus.pdf \\
        --perGroup

    plotHeatmap \\
        -m TSS_Plus_matrix.gz \\
        -o TSS_Metagene_Plus_Heatmap.pdf

    # Minus strand: computeMatrix and plot
    computeMatrix scale-regions \\
        -S \$MINUS_BAMS \\
        -R ${config.gtf} \\
        --beforeRegionStartLength ${window} \\
        --regionBodyLength 1000 \\
        --afterRegionStartLength ${window} \\
        --skipZeros \\
        -o TSS_Minus_matrix.gz \\
        -p 8

    plotProfile \\
        -m TSS_Minus_matrix.gz \\
        -o TSS_Metagene_Minus.pdf \\
        --perGroup

    plotHeatmap \\
        -m TSS_Minus_matrix.gz \\
        -o TSS_Metagene_Minus_Heatmap.pdf
    """
}
