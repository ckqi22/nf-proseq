process parse_config {

    output:
    path "config.txt"

    script:
    """
    python3 ${projectDir}/bin/parse_config.py \\
    --species_config ${params.species_config} \\
    --information_config ${params.information_config} \\
    --species ${params.species} \\
    ${params.build ? "--build ${params.build}" : ""} \\
    ${params.gtf ? "--gtf ${params.gtf}" : ""} \\
    ${params.genome_fasta ? "--genome_fasta ${params.genome_fasta}" : ""} \\
    --output config.txt
    """
}
