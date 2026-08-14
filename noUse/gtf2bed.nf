process GTF2BED {

    input:
    path gtf

    output:
    path '*.bed', emit: bed

    script:
    // This script is bundled with the pipeline, in bin/
    """
    gtf2bed \\
        ${gtf} \\
        > ${gtf.baseName}.bed
    """
}
