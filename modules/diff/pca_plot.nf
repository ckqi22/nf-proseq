process PCA_PLOT {

    input:
    path diff_dir

    output:
    path "pca_out/*", emit: pca

    script:
    """
    mkdir -p pca_out

    ${params.r} /workplace/pipeline/code/PCA_plot.R \\
      --input_dir ${diff_dir} \\
      --output_dir pca_out \\
      --width 7 \\
      --height 7 \\
      --dpi 600
    """
}
