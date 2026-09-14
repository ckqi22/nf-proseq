process SPIKEIN_SCALE {
    tag "spikein_scale"

    container "bio-base:1.0.0"

    // 收集全样本 *_spike_count.txt（每文件单行整数，样本名 = basename 去掉 .spike_count.txt），
    // 调 bin/spikein_scale.R 计算因子表。
    // 调用方式：SPIKEIN_SCALE(SPIKEIN_COUNT.out.counts.map{_m, f -> f}.collect())

    input:
    path spike_counts

    output:
    path "spikein_scale_factors.tsv", emit: factors   // sample spike_count factor size_factor

    script:
    """
    source /home/ck/miniconda3/bin/activate renv

    for f in ${spike_counts}; do
        s=\$(basename "\$f" .spike_count.txt)
        printf '%s\t%s\n' "\$s" "\$(cat "\$f")" >> spike_counts.tsv
    done

    Rscript ${projectDir}/bin/spikein_scale.R \\
        --counts spike_counts.tsv \\
        --output_dir ./
    """
}
