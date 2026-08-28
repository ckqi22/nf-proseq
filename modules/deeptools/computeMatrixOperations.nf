process COMPUTEMATRIXOPERATIONS {
    tag "${meta.sample}"

    input:
    tuple val(meta), path(plus_matrix), path(minus_matrix)

    output:
    tuple val(meta), path("${meta.sample}_TSS_merged_matrix.gz"), emit: matrix

    script:
    // Uses deepTools built-in computeMatrixOperations rbind to merge two
    // strand-specific computeMatrix output files by row-concatenation.
    // This correctly handles the @-prefixed JSON header format and updates
    // group_boundaries / group_labels. Requires deepTools >= 3.1.2
    // (v3.1.1 had a bug where rbind did not update group_labels).
    """
    source /workplace/hanguojun/mambaforge/bin/activate deeptools

    computeMatrixOperations rbind \\
        -m ${plus_matrix} ${minus_matrix} \\
        -o ${meta.sample}_TSS_merged_matrix.gz
    """
}
