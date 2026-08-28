process PROCESS_FASTP {
    tag "process_fastp"

    input:
    path json_files

    output:
    path "*.{tiff,pdf}", emit: base_quality_plot
    path "read_statistics.txt", emit: statistics

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
        --data_percentage ${params.threshold.data_percentage}

    ${params.r} /workplace/pipeline/code/base_quality_plot.R \\
        -i ./ \\
        -o ./ \\
        --width 7 --height 7 \\
        --colour "#C85D4D" \\
        --dpi 600
    """
}