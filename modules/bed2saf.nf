process BED2SAF {
    tag "${meta.id}"

    input:
    tuple val(meta), path(bed)

    output:
    tuple val(meta), path("*.saf"), emit: saf

    script:
    """
    awk 'OFS="\\t" {print \$1"."\$2"."\$3, \$1, \$2, \$3, "."}' \\
        ${bed} \\
        > ${bed.baseName}.saf
    """
}
