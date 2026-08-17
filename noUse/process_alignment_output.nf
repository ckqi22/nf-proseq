process process_alignment_output {
    tag "${meta.sample}"

    input:
    tuple val(meta), path(genomeRate), path(alignment_log)

    output:
    path "${meta.sample}.alignment_stats.txt", emit: stats

    script:
    """
    # Parse Bowtie2 summary and samtools log to extract alignment statistics
    echo "=== Alignment Statistics: ${meta.sample} ===" > ${meta.sample}.alignment_stats.txt
    echo "" >> ${meta.sample}.alignment_stats.txt

    # Bowtie2 overall alignment rate
    if [ -f "${genomeRate}" ]; then
        echo "--- Bowtie2 Summary ---" >> ${meta.sample}.alignment_stats.txt
        cat ${genomeRate} >> ${meta.sample}.alignment_stats.txt
        echo "" >> ${meta.sample}.alignment_stats.txt

        # Extract key metrics
        overall=\$(grep "overall alignment rate" ${genomeRate} | head -1)
        echo "Overall: \$overall" >> ${meta.sample}.alignment_stats.txt
    fi

    # Samtools stats from alignment log
    if [ -f "${alignment_log}" ]; then
        echo "" >> ${meta.sample}.alignment_stats.txt
        echo "--- Samtools Log ---" >> ${meta.sample}.alignment_stats.txt
        tail -5 ${alignment_log} >> ${meta.sample}.alignment_stats.txt
    fi
    """
}
