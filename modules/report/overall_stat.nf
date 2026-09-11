process OVERALL_STAT {
    tag "overall_stat"

    input:
    path read_statistics      // read_statistics.txt（PROCESS_FASTP 产物，单条）
    path bowtie2_summaries    // *_summary_bowtie2.txt（collect 后的 List，整体平铺进 workdir）

    output:
    path "overall_statistics.txt", emit: overall_statistics

    script:
    """
    ${params.r} ${projectDir}/bin/report/overall_statistics.R \\
        --read_statistics ${read_statistics} \\
        --bowtie2_summary_dir . \\
        --output overall_statistics.txt
    """
}
