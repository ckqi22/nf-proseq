#!/usr/bin/env nextflow
//
// ============================================================
// nf_PROseq — Nextflow DSL2 PRO-seq Analysis Pipeline
// ============================================================
//
// dUTP strand-specific library | --fr alignment | featureCounts -s 2
//
// Workflow:
//   1. Fastp QC + adapter trimming
//   2. Bowtie2 --fr alignment + strand split
//   3. featureCounts: full gene / TSS pause / gene body (3 regions per gene)
//   4. Pausing index (TSS / gene body)
//   5. Metagene TSS/TES profiles (deepTools)
//
// Samplesheet: sample,group,r1,r2
// ============================================================

include { parse_config      } from './modules/parse_config.nf'
include { preprocess        } from './subworkflows/preprocess.nf'
include { quantify          } from './subworkflows/quantify.nf'
include { pause_analysis    } from './subworkflows/pause_analysis.nf'
include { metagene_analysis } from './subworkflows/metagene_analysis.nf'

workflow {

    // ========================================================================
    // Step 0: Parse genome configuration
    // ========================================================================
    config_ch = parse_config().map { text ->
        def config = [:]
        text.split('\n').each { line ->
            def p = line.split(': ', 2)
            if (p.size() >= 2) { config[p[0].trim()] = p[1].trim() }
        }
        println "============================================"
        println " PRO-seq pipeline — config parsed"
        println " bowtie2_index : ${config.bowtie2_index}"
        println " gtf           : ${config.gtf}"
        println " build         : ${config.build}"
        println " sample_sheet  : ${params.sample_sheet}"
        println "============================================"
        return config
    }

    // ========================================================================
    // Step 1: Read samplesheet
    // ========================================================================
    read_ch = Channel.fromPath(params.sample_sheet)
        .splitCsv(header: true)
        .map { row ->
            if (!row.sample)        { error "samplesheet missing 'sample' column" }
            if (!row.r1 || !row.r2) { error "samplesheet missing 'r1' or 'r2' column" }
            def meta = [
                sample: row.sample,
                group:  row.group ?: 'unknown'
            ]
            [meta, [file(row.r1), file(row.r2)]]
        }

    // ========================================================================
    // Step 2: Preprocessing (Fastp + Bowtie2 --fr + strand split)
    // ========================================================================
    preprocess(read_ch, params.adapter ?: 'I', config_ch)

    // ========================================================================
    // Step 3: Quantification (full gene + TSS + gene body, single file)
    // ========================================================================
    quantify(preprocess.out.bam, config_ch)

    // ========================================================================
    // Step 4: Pausing index (TSS / gene body)
    // ========================================================================
    groups_config_ch = Channel.value(params.group ?: [:])
    pause_analysis(quantify.out.counts, groups_config_ch)

    // ========================================================================
    // Step 5: Metagene (TSS/TES profiles via deepTools)
    // ========================================================================
    plus_bam_list = preprocess.out.strand_bams
        .map { _meta, pbam, _mbam -> pbam }
        .collectFile(name: 'plus_bams.txt', newLine: true) { "${it}\n" }

    minus_bam_list = preprocess.out.strand_bams
        .map { _meta, _pbam, mbam -> mbam }
        .collectFile(name: 'minus_bams.txt', newLine: true) { "${it}\n" }

    metagene_analysis(plus_bam_list, minus_bam_list, config_ch)

    // ========================================================================
    // Publish results to output directories
    // ========================================================================
    publish:
    cutadapt_json      = preprocess.out.cutadapt_json
    cutadapt_html      = preprocess.out.cutadapt_html
    cutadapt_log       = preprocess.out.cutadapt_log
    trimmed_reads      = preprocess.out.trimmed_reads
    base_quality_plot  = preprocess.out.base_quality_plot
    statistics         = preprocess.out.statistics
    statistics_log     = preprocess.out.statistics_log

    bam                = preprocess.out.bam
    bai                = preprocess.out.bai
    genomeRate         = preprocess.out.genomeRate
    alignment_log      = preprocess.out.alignment_log

    counts             = quantify.out.counts
    summary            = quantify.out.summary

    pi_all             = pause_analysis.out.pi_all
    pi_boxplot         = pause_analysis.out.pi_boxplot
    pi_diff            = pause_analysis.out.pi_diff

    tss_plus_pdf       = metagene_analysis.out.tss_plus_pdf
    tss_minus_pdf      = metagene_analysis.out.tss_minus_pdf
    tes_plus_pdf       = metagene_analysis.out.tes_plus_pdf
    tes_minus_pdf      = metagene_analysis.out.tes_minus_pdf
}

// ------------------------------------------------------------------
// Output directive
// ------------------------------------------------------------------
output {
    cutadapt_json     { path "03.Data_QC/" }
    cutadapt_html     { path "03.Data_QC/" }
    cutadapt_log      { path "03.Data_QC/" }
    trimmed_reads     { path "03.Data_QC/" }
    base_quality_plot { path "03.Data_QC/" }
    statistics        { path "03.Data_QC/" }
    statistics_log    { path "03.Data_QC/" }

    bam               { path "04.Alignment/" }
    bai               { path "04.Alignment/" }
    genomeRate        { path "04.Alignment/" }
    alignment_log     { path "04.Alignment/" }

    counts            { path "05.Quantification/" }
    summary           { path "05.Quantification/" }

    pi_all            { path "06.Pausing_Index/" }
    pi_boxplot        { path "06.Pausing_Index/" }
    pi_diff           { path "06.Pausing_Index/" }

    tss_plus_pdf      { path "07.Metagene/" }
    tss_minus_pdf     { path "07.Metagene/" }
    tes_plus_pdf      { path "07.Metagene/" }
    tes_minus_pdf     { path "07.Metagene/" }
}
