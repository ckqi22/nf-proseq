process process_cutadapt_output {
    tag "${meta.sample}"

    input:
    tuple val(meta), path(json_files)   // 来自 cutadapt 进程的所有 JSON 文件列表

    output:
    path "*.{tiff,pdf}", emit: base_quality_plot
    path "read_statistics.txt", emit: statistics   // 最终汇总文件
    path "z.read_statistics.log", emit: statistics_log

    script:
    """
    mkdir -p stat
    
    source /workplace/hanguojun/mambaforge/bin/activate snakemake
    python3 /workplace/pipeline/WTSS/scripts/fastp_results_organise.py \\
        -i ./ \\
        -o ./ \\
        --q30_thres 0.7 \\
        --dup_thres ${params.threshold.duplicate_rate} \\
        --data_amount ${params.threshold.data_amount} \\
        --data_percentage ${params.threshold.data_percentage} \\
        > z.read_statistics.log 2>&1

    ${params.r} /workplace/pipeline/code/base_quality_plot.R \\
        -i ./ \\
        -o ./ \\
        --width 7 --height 7 \\
        --colour "#C85D4D" \\
        --dpi 600 \\
        >> z.read_statistics.log 2>&1
    """
}