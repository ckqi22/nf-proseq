process DE_PLOT {

    input:
    path diff_dir

    output:
    path "de_plot/*", emit: plot

    script:
    """
    mkdir -p de_plot

    ${params.r} /workplace/pipeline/code/scatter_plot.R \\
        --input_dir ${diff_dir} \\
        --output_dir de_plot \\
        --width 7 \\
        --height 7 \\
        --dpi 600

    ${params.r} /workplace/pipeline/code/volcano_plot.R \\
        --input_dir ${diff_dir} \\
        --output_dir de_plot \\
        --width 7 \\
        --height 5 \\
        --dpi 600

    ${params.r} /workplace/pipeline/code/heatmap_plot.R \\
        --input_dir ${diff_dir} \\
        --output_dir de_plot \\
        --width 6 \\
        --height 8 \\
        --dpi 600
    """
}
