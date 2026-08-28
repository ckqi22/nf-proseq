process POL2_COUNT_MERGE {
    tag "pol2"

    container "bio-base:1.0.0"

    input:
    path(counts)          // list of per-sample *.pol2.counts.txt (collected)

    output:
    path "pol2.matrix.txt", emit: matrix

    script:
    """
    source /home/ck/miniconda3/bin/activate renv

    Rscript ${projectDir}/bin/pol2_count_merge.R \\
        --inputs ${counts.join(',')} \\
        --output pol2.matrix.txt
    """
}
