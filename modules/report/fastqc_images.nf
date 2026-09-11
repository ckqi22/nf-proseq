process FASTQC_IMAGES {
    tag "fastqc_images"

    input:
    path raw_zips,     stageAs: 'raw/*'      // collect 后的 List：逐文件进 raw/ 子目录
    path trimmed_zips, stageAs: 'trimmed/*'

    output:
    path "fastqc_images", emit: images        // 目录（含 raw/ + trimmed/）

    script:
    """
    bash ${projectDir}/bin/report/extract_fastqc_images.sh raw trimmed fastqc_images
    """
}
