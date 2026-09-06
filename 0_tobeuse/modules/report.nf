process REPORT {

    input:
    path diff_dir         // deseq2_out/
    path enrich_dir       // gokegg_result/ + gsea_result/
    path metagene_dir     // metagene 结果目录
    path pausing_dir      // Pausing_Index 结果目录
    path pol2_signal      // pol2_signal_table.tsv
    path plot_dir         // de_plots/
    path config_yml       // params.yml
    path samplesheet_csv  // samplesheet.csv
    path statistics       // read_statistics.txt

    output:
    path "report/*", emit: report

    script:
    """
    mkdir -p report

    ${params.r} ${workflow.projectDir}/0_tobeuse/scripts/Report.R \\
      --output_dir report \\
      --config ${config_yml} \\
      --samplesheet ${samplesheet_csv} \\
      --diff_dir ${diff_dir} \\
      --enrich_dir ${enrich_dir} \\
      --metagene_dir ${metagene_dir} \\
      --pausing_dir ${pausing_dir} \\
      --pol2_signal ${pol2_signal} \\
      --plot_dir ${plot_dir} \\
      --statistics ${statistics} \\
      --resources_dir /workplace/pipeline/resources
    """
}
