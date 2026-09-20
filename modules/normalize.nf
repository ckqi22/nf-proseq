process NORMALIZE {
    tag "${name}"

    container "bio-base:1.0.0"

    input:
    path matrix          // count matrix (gene_id, length, samples)
    val methods          // comma-separated: cpm,fpkm[,rpkm]
    val name             // profile name (genebody | promoter | pol2)
    path spike_factors   // spikein_scale_factors.tsv；空 list 表示无 spike（不追加 .Spike 列）
    path total_mapped    // collected *.total_mapped.txt（每样本单行 read1 mapped count）

    output:
    path "${name}.normalized.txt", emit: normalized

    script:
    def spike_arg = spike_factors ? "--spike_factors ${spike_factors}" : ""
    """
    source /home/ck/miniconda3/bin/activate renv

    Rscript ${projectDir}/bin/normalize.R \\
        --input ${matrix} \\
        --methods ${methods} \\
        --total_mapped ${total_mapped.join(',')} \\
        ${spike_arg} \\
        --output ${name}.normalized.txt
    """
}
