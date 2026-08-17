process FEATURECOUNTS_MERGE {
    tag "${type}"

    container "bio-base:1.0.0"

    input:
    path(counts)          // list of per-sample *.featureCounts.txt (collected)
    val type              // 'tss' | 'genebody'

    output:
    path "${type}.matrix.txt", emit: matrix

    script:
    """
    source /home/ck/miniconda3/bin/activate renv

    Rscript ${projectDir}/bin/featurecounts_merge.R \\
        --inputs ${counts.join(',')} \\
        --output ${type}.matrix.txt
    """
}
