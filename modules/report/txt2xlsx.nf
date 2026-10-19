process TXT2XLSX {
    tag "txt2xlsx"

    input:
    path read_statistics      // read_statistics.txt
    path overall_statistics   // overall_statistics.txt（OVERALL_STAT 产物）
    val  diff_dir             // deseq2_out/ 绝对路径字符串（无差异分析时 = ''）
    path profile_file         // genebody.annotated.txt（normalize.R 产物，main_mapped 口径表达谱）

    output:
    path "report_xlsx", emit: xlsx_dir        // 目录

    script:
    def cg = params.diff?.compared_groups ?: []
    def fc = ''
    def p  = ''
    if (cg) {
        def fcs = cg.collect { it.split(',')[2].trim() }.unique()
        def ps  = cg.collect { it.split(',')[3].trim() }.unique()
        fc = (fcs.size() == 1) ? fcs[0] : ''
        p  = (ps.size()  == 1) ? ps[0]  : ''
    }
    def de_args = (fc && p) ? "--fc ${fc} --p ${p}" : ""
    """
    set -euo pipefail
    mkdir -p report_xlsx

    # ---- 1. read_statistics ----
    ${params.r} ${projectDir}/bin/report/txt2xlsx.R \\
        --input ${read_statistics} \\
        --output report_xlsx/read_statistics.xlsx \\
        --title "Read Statistics" \\
        --sort_by Sample \\
        --note_file ${projectDir}/templates/read_statistics.txt

    # ---- 2. overall_statistics ----
    ${params.r} ${projectDir}/bin/report/txt2xlsx.R \\
        --input ${overall_statistics} \\
        --output report_xlsx/overall_statistics.xlsx \\
        --title "Overall Statistics" \\
        --sort_by Sample \\
        --note_file ${projectDir}/templates/overall_statistics.txt

    # ---- 3. expression profiling (GeneBody 表达谱；取 normalize.R 的 main_mapped 口径产物)----
    if [ -f "${profile_file}" ]; then
        ${params.r} ${projectDir}/bin/report/txt2xlsx.R \\
            --input ${profile_file} \\
            --output report_xlsx/expression_profiling.xlsx \\
            --title "nascent mRNA Expression Profiling" \\
            --title_align left \\
            --note_file ${projectDir}/templates/expression_profiling.txt
    fi

    # ---- 4/5. 差异表 (_all / _de)----
    if [ -n "${diff_dir}" ] && ls "${diff_dir}"/*_all.txt >/dev/null 2>&1; then
        ${params.r} ${projectDir}/bin/report/txt2xlsx.R \\
            --diff_dir "${diff_dir}" \\
            --diff_suffix _all.txt \\
            --output report_xlsx/All_Comparisons_genes.xlsx \\
            --group_by Regulation \\
            --sheet_strip 2 \\
            --title "All Comparisons (DEG statistics)" \\
            --note_file ${projectDir}/templates/all_comparisons.txt \\
            --profile_file ${profile_file}
    fi

    if [ -n "${diff_dir}" ] && ls "${diff_dir}"/*_de.txt >/dev/null 2>&1; then
        ${params.r} ${projectDir}/bin/report/txt2xlsx.R \\
            --diff_dir "${diff_dir}" \\
            --diff_suffix _de.txt \\
            --output report_xlsx/Differentially_Expressed_genes.xlsx \\
            --group_by Regulation \\
            --sheet_strip 2 \\
            --title "Differentially Expressed Genes" \\
            --note_file ${projectDir}/templates/differentially_expressed.txt \\
            --profile_file ${profile_file} \\
            ${de_args}
    fi
    """
}
