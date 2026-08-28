process FEATURECOUNTS {
    tag "${meta.sample}"

    input:
    tuple val(meta), path(bams)
    path annotation
    val type

    output:
    path "*.featureCounts.txt", emit: counts

    script:
    // PRO-seq 定量永远吃 R1-only 单端 BAM（R2 是 5' 接头侧、无 Pol II 信号），
    // 故不加 -p；否则 featureCounts 单端模式会因 paired flag 报错。
    def is_gtf = annotation.name.endsWith("gtf")
    def anno_fmt_args = is_gtf ? "-F GTF" : "-F SAF"
    def type_args = is_gtf ? "-t gene" : ""
    def attr_args = is_gtf ? "-g gene_id" : "-g GeneID" 

    // PRO-seq 标准库 read1 反义(antisense) → 计数取反链(-s 2)。
    // strandedness 已从 params.yml 移除、待议定；暂默认 reverse，议定后用
    // params.strandedness 覆盖（'forward'→-s 1, 'reverse'→-s 2, 'unstranded'→-s 0）。
    def strandedness = 0
    def lib_strand = params.strandedness ?: 'reverse'
    if (lib_strand == 'forward') {
        strandedness = 1
    } else if (lib_strand == 'reverse') {
        strandedness = 2
    }

    """
    source /workplace/hanguojun/mambaforge/bin/activate snakemake

    ${params.feature_counts} \\
        -T 2 \\
        ${anno_fmt_args} \\
        -a ${annotation} \\
        ${type_args} \\
        ${attr_args} \\
        -s ${strandedness} \\
        -o ${meta.sample}_${type}.featureCounts.txt \\
        ${bams.join(' ')}
    """
}
