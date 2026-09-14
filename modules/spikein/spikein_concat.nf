process SPIKEIN_CONCAT {

    input:
    val main_fasta
    val spike_fasta

    output:
    path "combined.fa",  emit: combined_fasta
    path "spike.chroms", emit: spike_chroms

    script:
    """
    gunzip -cf ${main_fasta} > combined.fa
    gunzip -cf ${spike_fasta} \\
        | awk -v prefix="spike_" '/^>/{print ">" prefix substr(\$0,2); next} {print}' \\
        >> combined.fa
    
    gunzip -cf ${spike_fasta} \\
        | awk -v prefix="spike_" '/^>/{sub(/^>/,""); print prefix \$1}' \\
        > spike.chroms
    """
}
