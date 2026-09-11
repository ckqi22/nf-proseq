process PAUSING_XLSX {
    tag "pausing_xlsx"

    input:
    path pi_table          // Pausing_Index.tsv（pausing_index.R 产物）
    path promoter_bed      // promoter.bed（BED6 无表头，代表 transcript TSS 窗口）
    path genebody_bed      // genebody.bed（BED6 无表头，代表 transcript gene body）
    val  annotation        // gene_annotation 路径字符串（可为 ''，内容为空跳过）

    output:
    path "PROSeq_pausing.xlsx", emit: xlsx

    script:
    def ann_arg = (annotation && annotation.toString().trim()) ? "--annotation ${annotation}" : ""
    def upstream        = params.tss.upstream        ?: 50
    def downstream      = params.tss.downstream      ?: 300
    def genebody_offset = params.tss.genebody_offset ?: 301
    """
    ${params.r} ${projectDir}/bin/report/pausing_xlsx.R \\
        --pi ${pi_table} \\
        --promoter_bed ${promoter_bed} \\
        --genebody_bed ${genebody_bed} \\
        ${ann_arg} \\
        --tss_upstream ${upstream} \\
        --tss_downstream ${downstream} \\
        --genebody_offset ${genebody_offset} \\
        --output PROSeq_pausing.xlsx
    """
}
