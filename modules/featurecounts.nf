process FEATURECOUNTS {
    tag "${meta.sample}"

    input:
    tuple val(meta), path(bams)
    path annotation
    val type

    output:
    path "*.featureCounts.txt", emit: counts

    script:
    def paired_end_args = meta.single_end ? "" : "-p"
    def is_gtf = annotation.name.endsWith("gtf")
    def anno_fmt_args = is_gtf ? "-F GTF" : "-F SAF"
    def type_args = is_gtf ? "-t gene" : ""
    def attr_args = is_gtf ? "-g gene_id" : "-g GeneID" 

    def strandedness = 0
    if (params.strandedness == 'forward') {
        strandedness = 1
    } else if (params.strandedness == 'reverse') {
        strandedness = 2
    }

    """
    source /workplace/hanguojun/mambaforge/bin/activate snakemake

    ${params.feature_counts} \\
        ${paired_end_args} \\
        -T 2 \\
        ${anno_fmt_args} \\
        -a ${annotation} \\
        ${type_args} \\
        ${attr_args} \\
        -s 1 \\
        -o ${meta.sample}_${type}.featureCounts.txt \\
        ${bams.join(' ')}
    """
}
