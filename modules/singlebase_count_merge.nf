process SINGLEBASE_COUNT_MERGE {
    tag "${out_name}"

    container "bio-base:1.0.0"

    input:
    path(counts)          // list of per-sample *.${type}.counts.txt (collected)
    val type              // region type — used to derive the file suffix (e.g. promoter)
    val out_name          // output matrix base name (e.g. genebody_perbase)

    output:
    path "${out_name}.matrix.txt", emit: matrix

    script:
    """
    source /home/ck/miniconda3/bin/activate renv

    Rscript ${projectDir}/bin/singlebase_count_merge.R \\
        --inputs ${counts.join(',')} \\
        --suffix .${type}.counts.txt \\
        --output ${out_name}.matrix.txt
    """
}
