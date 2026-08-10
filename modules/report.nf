process report {
    tag "report"

    input:
    path all_results

    output:
    path "PRO-seq_Report.pdf", emit: report_pdf

    script:
    """
    source /workplace/hanguojun/mambaforge/bin/activate snakemake

    # Generate comprehensive PDF report using R Markdown

    ${params.r} -e "rmarkdown::render('${projectDir}/bin/proseq_report.Rmd', \\
        output_file='PRO-seq_Report.pdf', \\
        params=list(results_dir='${all_results}'))"
    """
}
