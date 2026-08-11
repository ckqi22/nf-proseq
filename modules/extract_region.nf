process extract_region {
    tag "${region}"

    input:
    path(count_files)
    val(region)       // _tss, _gene_body, _full_gene

    output:
    path "merged_${region}.txt", emit: matrix

    script:
    def flist = count_files.join(',')
    """
    ${params.r} ${projectDir}/bin/extract_region_matrix.R \\
        --files ${flist} \\
        --region ${region} \\
        --output merged_${region}.txt
    """
}
