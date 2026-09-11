process TXT2XLSX {
    tag "txt2xlsx"

    input:
    path read_statistics      // read_statistics.txt
    path overall_statistics   // overall_statistics.txt（OVERALL_STAT 产物）
    val  diff_dir             // deseq2_out/ 绝对路径字符串（无差异分析时 = ''）

    output:
    path "report_xlsx", emit: xlsx_dir        // 目录

    script:
    """
    set -euo pipefail
    mkdir -p report_xlsx

    # ---- 1. read_statistics ----
    ${params.r} ${projectDir}/bin/report/txt2xlsx.R \\
        --input ${read_statistics} \\
        --output report_xlsx/read_statistics.xlsx \\
        --title "Read Statistics" \\
        --sort_by Sample \\
        --note "Fastp 质控统计：Raw Reads / Clean Reads 为过滤前后 reads 数，Clean Ratio 为保留比例，Q30 为 Q30 碱基比例，Duplication Rate 为重复率。"

    # ---- 2. overall_statistics ----
    ${params.r} ${projectDir}/bin/report/txt2xlsx.R \\
        --input ${overall_statistics} \\
        --output report_xlsx/overall_statistics.xlsx \\
        --title "Overall Statistics" \\
        --sort_by Sample \\
        --note "综合统计：读段质控与比对率汇总。"

    # ---- 3. expression profiling（GeneBody 表达谱）----
    if [ -n "${diff_dir}" ] && [ -f "${diff_dir}/GeneBody_expression_profiling.txt" ]; then
        ${params.r} ${projectDir}/bin/report/txt2xlsx.R \\
            --input "${diff_dir}/GeneBody_expression_profiling.txt" \\
            --output report_xlsx/expression_profiling.xlsx \\
            --title "nascent mRNA Expression Profiling" \\
            --title_align left \\
            --note "Gene Body 新生 RNA 表达谱：.Count 为原始 read 计数，.Fpkm 为归一化表达量（Fragments Per Kilobase of gene per Million mapped reads）。"
    fi

    # ---- 4/5. 差异表（_all / _de）----
    if [ -n "${diff_dir}" ] && ls "${diff_dir}"/*_all.txt >/dev/null 2>&1; then
        ${params.r} ${projectDir}/bin/report/txt2xlsx.R \\
            --diff_dir "${diff_dir}" \\
            --diff_suffix _all.txt \\
            --output report_xlsx/All_Comparisons_genes.xlsx \\
            --group_by Regulation \\
            --sheet_strip 2 \\
            --title "All Comparisons (DEG statistics)" \\
            --note "差异分析：FoldChange 为组间倍数，log2FoldChange 为 log2 倍数，Pvalue/FDR 为显著性，Regulation 为 up/down。"
    fi

    if [ -n "${diff_dir}" ] && ls "${diff_dir}"/*_de.txt >/dev/null 2>&1; then
        ${params.r} ${projectDir}/bin/report/txt2xlsx.R \\
            --diff_dir "${diff_dir}" \\
            --diff_suffix _de.txt \\
            --output report_xlsx/Differentially_Expressed_genes.xlsx \\
            --group_by Regulation \\
            --sheet_strip 2 \\
            --title "Differentially Expressed Genes" \\
            --note "差异基因：满足阈值（|log2FoldChange| 与 FDR 由 config 的 compared_groups 决定）的显著差异基因。"
    fi
    """
}
