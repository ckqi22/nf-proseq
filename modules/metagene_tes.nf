process metagene_tes {
    tag "metagene_tes"

    input:
    path plus_bam_list
    path minus_bam_list
    val config

    output:
    path "TES_Metagene_Plus.pdf",   emit: tes_plus_pdf
    path "TES_Metagene_Minus.pdf",  emit: tes_minus_pdf

    script:
    def window = params.tes_metagene_window ?: 3000

    """
    source /workplace/hanguojun/mambaforge/bin/activate snakemake

    # TES metagene plot (strand-specific) using deepTools
    # Plus strand BAMs; Minus strand BAMs

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
        -o TES_Plus_matrix.gz \\
        -p 8

    plotProfile \\
        -m TES_Plus_matrix.gz \\
        -o TES_Metagene_Plus.pdf \\
        --perGroup

    plotHeatmap \\
        -m TES_Plus_matrix.gz \\
        -o TES_Metagene_Plus_Heatmap.pdf

    # Minus strand: computeMatrix and plot
    computeMatrix scale-regions \\
        -S \$MINUS_BAMS \\
        -R ${config.gtf} \\
        --beforeRegionStartLength ${window} \\
        --regionBodyLength 1000 \\
        --afterRegionStartLength ${window} \\
        --skipZeros \\
        -o TES_Minus_matrix.gz \\
        -p 8

    plotProfile \\
        -m TES_Minus_matrix.gz \\
        -o TES_Metagene_Minus.pdf \\
        --perGroup

    plotHeatmap \\
        -m TES_Minus_matrix.gz \\
        -o TES_Metagene_Minus_Heatmap.pdf
    """
}
