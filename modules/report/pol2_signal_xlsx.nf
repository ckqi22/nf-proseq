process POL2_SIGNAL_XLSX {
    tag "pol2_signal_xlsx"

    input:
    path tsv          // pol2_signal_table.tsv（signal_table.R 产物，无 note 头）
    path note_file    // pol2_signal_table.note.txt（signal_table.R 侧车 note）

    output:
    path "PROSeq_pol2_signal.xlsx", emit: xlsx

    script:
    """
    ${params.r} ${projectDir}/bin/report/txt2xlsx.R \\
        --input ${tsv} \\
        --output PROSeq_pol2_signal.xlsx \\
        --title "Pol II Active Site Signal" \\
        --note_file ${note_file}
    """
}
