process igv_tracks {
    tag "igv_tracks"

    input:
    path bam_files
    path bigwig_files

    output:
    path "IGV_tracks", emit: igv_output

    script:
    """
    source /workplace/hanguojun/mambaforge/bin/activate snakemake

    mkdir -p IGV_tracks

    # Copy all BAM and BigWig files to the IGV_tracks directory
    # Organize by sample for easy loading in IGV

    ${params.r} ${projectDir}/bin/igv_tracks.R \\
        --bam_dir . \\
        --bw_dir . \\
        --output_dir IGV_tracks
    """
}
