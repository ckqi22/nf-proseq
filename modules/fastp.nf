process FASTP {
    tag "${meta.sample}"

    input:
    tuple val(meta), path(reads)
    val adapter_type

    output:
    tuple val(meta), path("*_R{1,2}_trimmed.fastq.gz"), emit: trimmed_reads
    tuple val(meta), path("*.fastp.json"), emit: json
    path("*.fastp.html"), emit: html
    path("*.log"), emit: log

    script:
    def args = ""
    if (adapter_type == "UMI" || adapter_type == "HT" || adapter_type == "SP") {
        args = "-U --umi_loc=per_read --umi_len=${params.umi_length} --umi_prefix=UMI --trim_poly_x"
    } else if (adapter_type == "I") {
        args = "" 
    }
    if (meta.single_end) {
        reads_args          = "--in1 ${reads}"
        trimmed_reads_args  = "--out1 ${meta.sample}_R1_trimmed.fastq.gz"
        detect_adapter_args = ""
    } else {
        reads_args          = "--in1 ${reads[0]} --in2 ${reads[1]}"
        trimmed_reads_args  = "--out1 ${meta.sample}_R1_trimmed.fastq.gz --out2 ${meta.sample}_R2_trimmed.fastq.gz"
        detect_adapter_args = "--detect_adapter_for_pe"
    }
    """
    source /workplace/hanguojun/mambaforge/bin/activate snakemake
    
    fastp \\
    ${reads_args} \\
    ${trimmed_reads_args} \\
    --json ${meta.sample}.fastp.json \\
    --html ${meta.sample}.fastp.html \\
    --report_title ${meta.sample} \\
    ${detect_adapter_args} \\
    ${args} > ${meta.sample}.log 2>&1
    """
}