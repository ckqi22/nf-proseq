#!/usr/bin/env nextflow
//
// ============================================================
// PRO-seq Nextflow Pipeline — Main Workflow
// ============================================================
// 14-step automated analysis for PRO-seq (Precision Run-On sequencing)
//
// Plan Reference:
//  01. Manual: project info
//  02. Manual: QC info
//  03. fastp QC -> bowtie2 alignment -> strand separation
//  04. featureCounts (gene body + TSS) -> spike correction -> strand-specific quantification
//  05. TSS pausing index calculation
//  06. edgeR differential expression
//  07. GO enrichment
//  08. KEGG enrichment
//  09. GSEA
//  10. Heatmap
//  11. Scatter plot
//  12. Volcano plot
//  13. Pol II active site profiling + metagene (TSS/TES)
//  14. bigwig + IGV + GEO + PDF report
//
// Samplesheet format: sample,group,strand,spike_ratio,R1,R2
// ============================================================

// ------------------------------------------------------------------
// Import subworkflows (NOT direct modules)
// ------------------------------------------------------------------
include { preprocess       } from './subworkflows/preprocess.nf'
include { quantify         } from './subworkflows/quantify.nf'
include { pause_analysis   } from './subworkflows/pause_analysis.nf'
include { diff_analysis    } from './subworkflows/diff_analysis.nf'
include { metagene_analysis } from './subworkflows/metagene_analysis.nf'
include { deliver          } from './subworkflows/deliver.nf'

// Import parse_config directly (needed for multi-subworkflow shared config)
include { parse_config     } from './modules/parse_config.nf'

// ------------------------------------------------------------------
// Main workflow
// ------------------------------------------------------------------
workflow {

    main:
    // ==================================================================
    // Step 0: Parse configuration from params.yml via parse_config
    // ==================================================================
    config_ch = parse_config().map { text ->
        def config = [:]
        text.split('\n').each { line ->
            def p = line.split(': ', 2)
            if (p.size() >= 2) {
                config[p[0]] = p[1]
            }
        }
        println "============================================"
        println " PRO-seq pipeline — config parsed"
        println " bowtie2_index : ${config.bowtie2_index}"
        println " gtf           : ${config.gtf}"
        println " gene_annotation: ${config.gene_annotation}"
        println " rRNA_index    : ${config.rRNA_index}"
        println " build         : ${config.build}"
        println " sample_sheet  : ${params.sample_sheet}"
        println " adapter       : ${params.adapter}"
        println "============================================"
        return config
    }

    // ==================================================================
    // Step 1: Read samplesheet
    // Format: sample, group, strand, spike_ratio, R1, R2
    // ==================================================================
    read_ch = channel.fromPath(params.sample_sheet)
                .splitCsv(header: true)
                .map { row ->
                    // Validate required columns
                    if (!row.sample) { error "samplesheet missing 'sample' column" }
                    if (!row.r1 || !row.r2) { error "samplesheet missing 'r1' or 'r2' column" }

                    def meta = [
                        sample: row.sample,
                        group:  row.group ?: 'unknown',
                        strand: row.strand ?: 'unstranded',
                        spike_ratio: row.spike_ratio ?: '1'
                    ]
                    println "[samplesheet] sample=${meta.sample}  group=${meta.group}  strand=${meta.strand}  spike_ratio=${meta.spike_ratio}"
                    [meta, [file(row.r1), file(row.r2)]]
                }

    // ==================================================================
    // Step 2: Preprocessing (03. Data QC + Alignment + Strand split)
    // ==================================================================
    preprocess(read_ch, params.adapter ?: 'I', config_ch)

    // ==================================================================
    // Step 3: Quantification (04. Gene body + TSS + Spike correction)
    // ==================================================================
    quantify(preprocess.out.strand_bams, config_ch)

    // ==================================================================
    // Step 4: Pausing Index (05. TSS pausing index calculation)
    // ==================================================================
    // Prepare groups configuration from params.yml
    // groups_config is a map: group_name -> [sample1, sample2, ...]
    groups_config_ch = channel.value(params.group ?: [:])

    pause_analysis(
        quantify.out.tss_counts,
        quantify.out.gene_body_counts,
        groups_config_ch
    )

    // ==================================================================
    // Step 5: Differential analysis (06-12. DE + Enrichment + Plots)
    // ==================================================================
    // Collect gene body counts into a merged matrix
    // Fork channel: used by both diff_analysis and deliver
    quantify.out.gene_body_counts
        .flatten()
        .into { gb_matrix_de; gb_matrix_deliver }

    // Build comparisons string from params compared-groups
    comparisons_ch = channel.value(params.compared_groups ?: [])

    diff_analysis(
        gb_matrix_de,
        comparisons_ch,
        config_ch
    )

    // ==================================================================
    // Step 6: Metagene analysis (13. TSS/TES metagene profiles)
    // ==================================================================
    // Prepare plus and minus BAM lists for metagene
    // strand_bams contains: tuple val(meta), path(plus_bam), path(minus_bam)
    plus_bam_list = preprocess.out.strand_bams
        .map { _meta, pbam, _mbam -> pbam }
        .collectFile(
            name: 'plus_bams.txt',
            newLine: true
        ) { bam -> "${bam}\n" }

    minus_bam_list = preprocess.out.strand_bams
        .map { _meta, _pbam, mbam -> mbam }
        .collectFile(
            name: 'minus_bams.txt',
            newLine: true
        ) { bam -> "${bam}\n" }

    metagene_analysis(plus_bam_list, minus_bam_list, config_ch)

    // ==================================================================
    // Step 7: Deliver (13-14. Pol II profiling + BigWig + IGV + GEO + Report)
    // ==================================================================
    deliver(
        preprocess.out.strand_bams,
        preprocess.out.bam,
        gb_matrix_deliver,
        channel.value(file(params.sample_sheet)),
        config_ch
    )

    // ==================================================================
    // Results are published via the output {} directive below
    // ==================================================================
}

// ------------------------------------------------------------------
// Output directive: publish results to organized directories
// ------------------------------------------------------------------
output {
    '03.Data_QC/cutadapt_json'     { path "03.Data_QC/" }
    '03.Data_QC/cutadapt_html'     { path "03.Data_QC/" }
    '03.Data_QC/cutadapt_log'      { path "03.Data_QC/" }
    '03.Data_QC/trimmed_reads'     { path "03.Data_QC/" }
    '03.Data_QC/base_quality_plot' { path "03.Data_QC/" }
    '03.Data_QC/statistics'        { path "03.Data_QC/" }
    '03.Data_QC/statistics_log'    { path "03.Data_QC/" }

    '04.Alignment/bam'           { path "04.Alignment/" }
    '04.Alignment/bai'           { path "04.Alignment/" }
    '04.Alignment/genomeRate'    { path "04.Alignment/" }
    '04.Alignment/alignment_log' { path "04.Alignment/" }

    '05.Quantification/GeneBody/plus_counts'  { path "05.Quantification/GeneBody/" }
    '05.Quantification/GeneBody/minus_counts' { path "05.Quantification/GeneBody/" }
    '05.Quantification/GeneBody/combined'     { path "05.Quantification/GeneBody/" }
    '05.Quantification/TSS/plus_counts'       { path "05.Quantification/TSS/" }
    '05.Quantification/TSS/minus_counts'      { path "05.Quantification/TSS/" }
    '05.Quantification/TSS/combined'          { path "05.Quantification/TSS/" }

    '06.Pausing_Index/pi_all'    { path "06.Pausing_Index/" }
    '06.Pausing_Index/pi_boxplot' { path "06.Pausing_Index/" }
    '06.Pausing_Index/pi_diff'   { path "06.Pausing_Index/" }

    '07.Diff/diff_genes'      { path "07.Diff/" }
    '07.Diff/all_comparisons' { path "07.Diff/" }

    '08.GO_KEGG/go_results'   { path "08.GO_KEGG/" }
    '08.GO_KEGG/kegg_results' { path "08.GO_KEGG/" }

    '09.GSEA/gsea_results'    { path "09.GSEA/" }

    '10.Plot/heatmap'  { path "10.Plot/" }
    '10.Plot/scatter'  { path "10.Plot/" }
    '10.Plot/volcano'  { path "10.Plot/" }

    '11.Metagene/tss_plus'  { path "11.Metagene/" }
    '11.Metagene/tss_minus' { path "11.Metagene/" }
    '11.Metagene/tes_plus'  { path "11.Metagene/" }
    '11.Metagene/tes_minus' { path "11.Metagene/" }

    '12.Pol_II_active_site/profiling'  { path "12.Pol_II_active_site/" }
    '12.Pol_II_active_site/tss_heatmap' { path "12.Pol_II_active_site/" }
    '12.Pol_II_active_site/bigwigs'    { path "12.Pol_II_active_site/" }
    '12.Pol_II_active_site/igv'        { path "12.Pol_II_active_site/" }
    '12.Pol_II_active_site/geo'        { path "12.Pol_II_active_site/" }

    'report' { path "report/" }
}
