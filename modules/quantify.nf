process quantify {
    
    input:
    path bam
    val config

    output:
    path "All.Sample.featureCounts.txt", emit: txt
    path "All.Sample.Transcript.featureCounts.txt", emit: transcript_txt
    path "All.Sample.Fragments.featureCounts.txt", emit: fragment_txt
    path "All.Sample.featureCounts.txt.summary", emit: summary
    path "All.Sample.Transcript.featureCounts.txt.summary", emit: transcript_summary
    path "All.Sample.Fragments.featureCounts.txt.summary", emit: fragment_summary
    path "*.bam.featureCounts.bam", emit: output_files
    path "*.{tiff,pdf}", emit: plot
    path "z.featureCount.log", emit: log

    script:
    def bam_files = bam.join(' ')
    """
    ${params.feature_counts} -p                                     -a ${config.gtf} -o All.Sample.featureCounts.txt            -T 4 ${bam_files} >           z.featureCount.log 2>&1
    ${params.feature_counts} -p --countReadPairs                    -a ${config.gtf} -o All.Sample.Fragments.featureCounts.txt  -T 4 ${bam_files} >>          z.featureCount.log 2>&1
    ${params.feature_counts} -p --countReadPairs -g "transcript_id" -a ${config.gtf} -o All.Sample.Transcript.featureCounts.txt -T 4 ${bam_files} >>          z.featureCount.log 2>&1
    ${params.feature_counts} -p                  -g "gene_id"       -a ${config.gtf} -o All.Sample.featureCounts.txt            -T 4 ${bam_files} -R BAM >>   z.featureCount.log 2>&1
    
    ${params.r} /workplace/pipeline/code/PCA_plot.R --featurecount All.Sample.featureCounts.txt --config ${launchDir}/params.yml -o ./ >> z.featureCount.log 2>&1    
    """
    
}
