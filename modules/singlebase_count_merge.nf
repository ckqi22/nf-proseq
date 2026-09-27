process SINGLEBASE_COUNT_MERGE {
    tag "${out_name}"

    container "bio-base:1.0.0"

    input:
    path(counts)          // list of per-sample count files (collected)
    val samples           // list of sample names (collected, same order as counts)
    val out_name          // output matrix base name (e.g. pol2_promoter)
    path tx2gene          // transcript_id -> gene_id map (gtf2bed.R tx2gene.tsv)

    output:
    path "${out_name}.matrix.txt", emit: matrix

    script:
    """
    source /home/ck/miniconda3/bin/activate renv

    Rscript ${projectDir}/bin/singlebase_count_merge.R \\
        --sample ${samples.join(',')} \\
        --count  ${counts.join(',')} \\
        --tx2gene ${tx2gene} \\
        --output ${out_name}.matrix.txt
    """
}
