process PLOT_HEATMAP {
    tag "${meta.sample}"

    container "bio-base:1.0.0"

    input:
    tuple val(meta), path(plus_bw), path(minus_bw)
    path tss_bed

    output:
    path("${meta.sample}_metagene_heatmap.png"), emit: plot
    path("${meta.sample}_metagene_heatmap_matrix.tsv"), emit: matrix

    script:
    def upstream   = params.metagene.heatmap_upstream   ?: 50
    def downstream = params.metagene.heatmap_downstream ?: 150
    """
    source /home/ck/miniconda3/bin/activate renv

    Rscript ${projectDir}/bin/metagene_heatmap.R \\
      --prefix ${meta.sample} \\
      --plus_bw  ${plus_bw} \\
      --minus_bw ${minus_bw} \\
      --tss_bed ${tss_bed} \\
      --upstream ${upstream} \\
      --downstream ${downstream} \\
      --bin 10 \\
      --outdir ./
    """
}
