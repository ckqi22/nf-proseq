process FEATURECOUNTS {
    tag "${meta.sample}"

    input:
    tuple val(meta), path(bams)
    tuple val(type), path(annotation)

    output:
    path "*.featureCounts.txt", emit: counts

    script:
    // PRO-seq 定量永远吃 R1-only 单端 BAM（R2 是 5' 接头侧、无 Pol II 信号），
    // 故不加 -p；否则 featureCounts 单端模式会因 paired flag 报错。
    def is_gtf = annotation.name.endsWith("gtf")
    def anno_fmt_args = is_gtf ? "-F GTF" : "-F SAF"
    def type_args = is_gtf ? "-t gene" : ""
    def attr_args = is_gtf ? "-g gene_id" : "-g GeneID" 

    // 链向由 params.strandedness 决定（与 genomecov 信号轨同一约定）：
    //   reverse（默认，标准 PRO-seq，R1 antisense）→ -s 2（反链计数）
    //   forward（R1 sense）→ -s 1（正链计数）
    // 不支持 unstranded（PRO-seq 信号必须有链向），非法值报错。
    def strand = params.strandedness?.trim() ?: 'reverse'
    if (strand != 'forward' && strand != 'reverse') {
        error "params.strandedness must be 'reverse' or 'forward', got: ${strand}" 
    }
    def strandedness = (strand == 'forward') ? 1 : 2
    """
    ${params.feature_counts} \\
        -T 2 \\
        ${anno_fmt_args} \\
        -a ${annotation} \\
        ${type_args} \\
        ${attr_args} \\
        -s ${strandedness} \\
        -O \\
        -o ${meta.sample}_${type}.featureCounts.txt \\
        ${bams.join(' ')}
    """
}
