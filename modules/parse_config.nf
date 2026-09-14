process parse_config {

    output:
    path "config.txt", emit: info

    script:
    """
    python3 ${projectDir}/bin/parse_config.py \\
    --species_config ${params.species_config} \\
    --information_config ${params.information_config} \\
    --species ${params.species} \\
    ${params.build ? "--build ${params.build}" : ""} \\
    ${params.gtf ? "--gtf ${params.gtf}" : ""} \\
    ${params.genome_fasta ? "--genome_fasta ${params.genome_fasta}" : ""} \\
    ${params.spike_genome ? "--spike_genome ${params.spike_genome}" : ""} \\
    ${params.spike_fasta ? "--spike_fasta ${params.spike_fasta}" : ""} \\
    ${params.spike_gtf ? "--spike_gtf ${params.spike_gtf}" : ""} \\
    ${params.spike_index ? "--spike_index ${params.spike_index}" : ""} \\
    ${params.spike_chroms ? "--spike_chroms ${params.spike_chroms}" : ""} \\
    --output config.txt
    """
}
