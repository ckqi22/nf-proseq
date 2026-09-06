process PLOT_PROFILE {
    tag "${meta.sample}"

    container "bio-base:1.0.0"

    input:
    tuple val(meta), path(plus_bw), path(minus_bw)
    path tss_bed

    output:
    path("${meta.sample}_metagene_profile.png"), emit: plot
    path("${meta.sample}_metagene_profile_matrix.tsv"), emit: matrix

    script:
    """
    source /home/ck/miniconda3/bin/activate renv

    Rscript ${projectDir}/bin/metagene_profile.R \\
      --mirror \\
      --prefix ${meta.sample} \\
      --plus_bw  ${plus_bw} \\
      --minus_bw ${minus_bw} \\
      --tss_bed ${tss_bed} \\
      --window ${params.metagene.tss_window} \\
      --bin 10 \\
      --outdir ./
    """
}