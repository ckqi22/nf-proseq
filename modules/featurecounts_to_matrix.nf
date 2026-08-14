process FEATURECOUNTS_TO_MATRIX {

    input:
    path counts          // featureCounts output (Geneid, Chr, Start, End, Strand, Length, <sample>.bam ...)

    output:
    path "matrix.txt", emit: matrix

    script:
    """
    source /workplace/hanguojun/mambaforge/bin/activate snakemake

    ${params.r} ${projectDir}/bin/featurecounts_to_matrix.R \\
        --input ${counts} \\
        --output matrix.txt
    """
}
