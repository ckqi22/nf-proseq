process pol2_profiling {
    tag "${meta.sample}"

    input:
    tuple val(meta), path(bam)
    path tss_bed           // Optional: TSS regions BED for heatmap

    output:
    path "${meta.sample}.PROSeq_profiling.xlsx",  emit: profiling_xlsx
    path "${meta.sample}.TSS_Heatmap.pdf",         emit: tss_heatmap

    script:
    def tss_arg = (tss_bed && tss_bed.name != 'null') ? "--tss_bed ${tss_bed}" : ""

    """
    source /workplace/hanguojun/mambaforge/bin/activate snakemake

    mkdir -p bam_dir_${meta.sample}
    ln -sf \$(readlink -f ${bam}) bam_dir_${meta.sample}/${meta.sample}.bam
    if [ -f "${bam}.bai" ]; then
        ln -sf \$(readlink -f ${bam}.bai) bam_dir_${meta.sample}/${meta.sample}.bam.bai
    fi

    ${params.r} ${projectDir}/bin/pol2_profiling.R \\
        --bam_dir bam_dir_${meta.sample} \\
        ${tss_arg} \\
        --output_dir ./

    if [ -f PROSeq_profiling.xlsx ]; then
        mv PROSeq_profiling.xlsx ${meta.sample}.PROSeq_profiling.xlsx
    fi
    if [ -f PROSeq_TSS_Heatmap.pdf ]; then
        mv PROSeq_TSS_Heatmap.pdf ${meta.sample}.TSS_Heatmap.pdf
    fi
    """
}
