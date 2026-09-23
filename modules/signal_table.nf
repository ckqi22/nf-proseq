process SIGNAL_TABLE {
    tag "signal_table"

    container "bio-base:1.0.0"

    input:
    val manifest          // 6 列 TSV（无表头）：group \t sample \t plus_raw \t minus_raw \t plus_cpm \t minus_cpm（bigWig basename）
    path bigwig_files     // 全部 raw/cpm 正负链 bigWig 扁平列表（staged，供 R 按 basename 读取）
    path promoter_bed     // pause 窗口 BED6（gtf2bed.R promoter.bed，0-based 半开）
    path promoter_count   // pol2_promoter.matrix.txt（gene_id, length, samples...）
    path genebody_count   // pol2_genebody.matrix.txt
    path rep_gtf          // 代表转录本 GTF（longest_tx.gtf）→ transcriptid 列
    path annotation       // gene 注释表（首列 gene_id，其余列透传）
    val orientation       // 'reverse' | 'forward' 

    output:
    path "pol2_signal_table.tsv", emit: signal_table
    path "pol2_signal_table.note.txt", emit: note

    script:
    def st             = params.signal_table ?: [:]
    def active_frac    = st.active_frac      ?: 0.005
    def min_genebody   = params.tss.min_genebody_length ?: 800
    def peak_frac      = st.peak_frac        ?: 0.1
    def noise_quantile = st.noise_quantile   ?: 0.9
    def min_reps       = st.min_reps         ?: 2
    """
    source /home/ck/miniconda3/bin/activate renv

    printf '%s\n' '${manifest}' > manifest.tsv

    Rscript ${projectDir}/bin/signal_table.R \\
        --manifest manifest.tsv \\
        --promoter_bed ${promoter_bed} \\
        --promoter_count ${promoter_count} \\
        --genebody_count ${genebody_count} \\
        --rep_gtf ${rep_gtf} \\
        --annotation ${annotation} \\
        --active_frac ${active_frac} \\
        --min_genebody_length ${min_genebody} \\
        --peak_frac ${peak_frac} \\
        --noise_quantile ${noise_quantile} \\
        --min_reps ${min_reps} \\
        --orientation ${orientation} \\
        --output .
    """
}
