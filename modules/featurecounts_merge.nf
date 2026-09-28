process FEATURECOUNTS_MERGE {
    tag "${type}"

    container "nf-proseq:1.0.0"

    input:
    path(counts)          // list of per-sample *.featureCounts.txt (collected)
    val type              // promoter | genebody

    output:
    path "${type}.matrix.txt", emit: matrix

    script:
    """
    Rscript ${projectDir}/bin/featurecounts_merge.R \\
        --inputs ${counts.join(',')} \\
        --output ${type}.matrix.txt
    """
}
