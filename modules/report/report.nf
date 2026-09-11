process REPORT {

    input:
    val  enrich_files        // List<Path>：gokegg/gsea 子目录（无差异分析时 = 空 List）
    val  metagene_files      // List<Path>：metagene plot + matrix（必非空）
    val  pausing_files       // List<Path>：PROSeq_pausing.xlsx（单文件，必非空）
    val  plot_files          // List<Path>：heatmap/scatter/volcano（无差异分析时 = 空 List）
    val  base_quality_files  // List<Path>：*.tiff/*.pdf（过滤前/后碱基质量图，必非空）
    val  group_heatmap_files // List<Path>：组级 metagene heatmap（无分组时 = 空 List）
    path xlsx_dir            // report_xlsx/（TXT2XLSX 产物目录）
    path fastqc_images       // fastqc_images/（含 raw/ trimmed/）
    path config_yml          // params.yml
    path samplesheet_csv     // samplesheet.csv
    path resolved_config     // config.txt（parse_config 产物，补 params.yml 的 build/gtf）

    output:
    path "report/*", emit: report   // {name}_{projectNo}_{Species}_PRO-seq_Sequencing_Report_{date}/

    script:
    def enrich_cp_files = enrich_files.findAll { !it.toString().endsWith('GOOD_LUCK.txt') }
    def enrich_cp   = enrich_cp_files ? "cp -r ${enrich_cp_files.join(' ')} stage/enrich/" : "true"
    def metagene_cp = metagene_files ? "cp -r ${metagene_files.join(' ')} stage/metagene/" : "true"
    def pausing_cp  = pausing_files  ? "cp -r ${pausing_files.join(' ')} stage/pausing/"  : "true"
    def plot_cp     = plot_files     ? "cp -r ${plot_files.join(' ')} stage/plot/"        : "true"
    def baseq_cp    = base_quality_files ? "cp -r ${base_quality_files.join(' ')} stage/base_quality/" : "true"
    def group_hm_cp = group_heatmap_files ? "cp -r ${group_heatmap_files.join(' ')} stage/group_heatmaps/" : "true"
    """
    mkdir -p report stage/enrich stage/metagene stage/pausing stage/plot stage/base_quality stage/group_heatmaps

    ${enrich_cp}
    ${metagene_cp}
    ${pausing_cp}
    ${plot_cp}
    ${baseq_cp}
    ${group_hm_cp}

    ${params.r} ${projectDir}/bin/report/Report.R \\
      --output_dir report \\
      --config ${config_yml} \\
      --samplesheet ${samplesheet_csv} \\
      --xlsx_dir ${xlsx_dir} \\
      --fastqc_images ${fastqc_images} \\
      --base_quality_plot stage/base_quality \\
      --enrich_dir stage/enrich \\
      --metagene_dir stage/metagene \\
      --pausing_dir stage/pausing \\
      --plot_dir stage/plot \\
      --group_heatmaps stage/group_heatmaps \\
      --resolved_config ${resolved_config} \\
      --resources_dir ${params.report_resources_dir}
    """
}
