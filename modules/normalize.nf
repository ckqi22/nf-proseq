process NORMALIZE {
    tag "${name}"

    container "bio-base:1.0.0"

    input:
    path matrix          // count matrix (gene_id, length, samples)
    val methods          // comma-separated: cpm,fpkm[,rpkm]
    val name             // profile name (genebody | promoter | pol2)
    val spike_factors    // spikein_scale_factors.tsv 的绝对路径；'' 表示无 spike（走库大小）

    output:
    path "${name}.normalized.txt", emit: normalized

    script:
    def spike_arg = spike_factors ? "--spike_factors ${spike_factors}" : ""
    """
    source /home/ck/miniconda3/bin/activate renv

    Rscript ${projectDir}/bin/normalize.R \\
        --input ${matrix} \\
        --methods ${methods} \\
        ${spike_arg} \\
        --output ${name}.normalized.txt
    """
}
