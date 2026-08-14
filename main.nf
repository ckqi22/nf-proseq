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
include { align_bowtie2     } from './subworkflows/align_bowtie2.nf'

include { prepare_genome    } from './subworkflows/prepare_genome.nf'
include { quantification    } from './subworkflows/quantification.nf'
include { pause_analysis    } from './subworkflows/pause_analysis.nf'
include { metagene_analysis } from './subworkflows/metagene_analysis.nf'

workflow {

    main:
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
        println " genome_fasta  : ${config.genome_fasta}"
        println " bowtie2_index : ${config.bowtie2_index}"
        println " gtf           : ${config.gtf}"
        println " build         : ${config.build}"
        println " sample_sheet  : ${params.sample_sheet}"
        println "============================================"
        return config
    }

    // Resolve reference genome (fasta + bowtie2 index) and prepare
    // TSS + gene body SAF annotations from the reference GTF
    prepare_genome(config_ch)

    // ========================================================================
    // Step 1: Read samplesheet
    // ========================================================================
    read_ch = Channel.fromPath(params.sample_sheet)
        .splitCsv(header: true)
        .map { row ->
            if (!row.sample)    { error "samplesheet missing 'sample' column" }
            if (!row.r1)        { error "samplesheet missing 'r1' column" }
            def meta = [
                sample:     row.sample,
                group:      row.group ?: 'unknown',
                single_end: row.r2 ? false : true 
            ]
            [meta, [file(row.r1), file(row.r2)]]
        }

    // ========================================================================
    // Step 2: preprocess
    // ========================================================================
    preprocess(read_ch, params.adapter ?: 'I')

    // ========================================================================
    // Step 3: alignment
    // ========================================================================
    align_bowtie2(preprocess.out.trimmed_reads, prepare_genome.out.index, prepare_genome.out.fasta)
        
    // ========================================================================
    // Step 4: quantification (full gene + TSS + gene body, single file)
    // ========================================================================
    bam_list = align_bowtie2.out.bam
        .map { meta, bam -> [meta, bam] }
        .collect()
        .map { items -> [items[0][0], items.collect { it[1] }] }

    quantification(bam_list, prepare_genome.out.tss_saf, prepare_genome.out.genebody_saf)

    // ========================================================================
    // Step 5: Pausing index (TSS / gene body)
    // ========================================================================
    groups_config_ch = Channel.value(params.group ?: [:]).view()
    // pause_analysis(quantification.out.tss_counts, quantification.out.genebody_counts, groups_config_ch)

    // ========================================================================
    // Step 6: Metagene (TSS/TES profiles via deepTools)
    // ========================================================================
    // plus_bam_list = preprocess.out.strand_bams
    //     .map { _meta, pbam, _mbam -> pbam }
    //     .collectFile(name: 'plus_bams.txt', newLine: true) { "${it}\n" }

    // minus_bam_list = preprocess.out.strand_bams
    //     .map { _meta, _pbam, mbam -> mbam }
    //     .collectFile(name: 'minus_bams.txt', newLine: true) { "${it}\n" }

    // metagene_analysis(plus_bam_list, minus_bam_list, config_ch)

    // ========================================================================
    // Publish results to output directories
    // ========================================================================
    publish:
    fastp_json         = preprocess.out.fastp_json
    fastp_html         = preprocess.out.fastp_html
    fastp_log          = preprocess.out.fastp_log
    trimmed_reads      = preprocess.out.trimmed_reads
    base_quality_plot  = preprocess.out.base_quality_plot
    statistics         = preprocess.out.statistics
    statistics_log     = preprocess.out.statistics_log

    bam                = align_bowtie2.out.bam
    bai                = align_bowtie2.out.bai
    genomeRate         = align_bowtie2.out.genomeRate
    alignment_log      = align_bowtie2.out.alignment_log

    tss_counts         = quantification.out.tss_counts
    genebody_counts    = quantification.out.genebody_counts

    // pi_all             = pause_analysis.out.pi_all
    // pi_boxplot         = pause_analysis.out.pi_boxplot
    // pi_diff            = pause_analysis.out.pi_diff

    // tss_plus_pdf       = metagene_analysis.out.tss_plus_pdf
    // tss_minus_pdf      = metagene_analysis.out.tss_minus_pdf
    // tes_plus_pdf       = metagene_analysis.out.tes_plus_pdf
    // tes_minus_pdf      = metagene_analysis.out.tes_minus_pdf
}

// ------------------------------------------------------------------
// Output directive
// ------------------------------------------------------------------
output {
    fastp_json        { path "03.Data_QC/" }
    fastp_html        { path "03.Data_QC/" }
    fastp_log         { path "03.Data_QC/" }
    trimmed_reads     { path "03.Data_QC/" }
    base_quality_plot { path "03.Data_QC/" }
    statistics        { path "03.Data_QC/" }
    statistics_log    { path "03.Data_QC/" }

    bam               { path "04.Alignment/" }
    bai               { path "04.Alignment/" }
    genomeRate        { path "04.Alignment/" }
    alignment_log     { path "04.Alignment/" }

    tss_counts        { path "05.Quantification/" }
    genebody_counts   { path "05.Quantification/" }

    // pi_all            { path "06.Pausing_Index/" }
    // pi_boxplot        { path "06.Pausing_Index/" }
    // pi_diff           { path "06.Pausing_Index/" }

    // tss_plus_pdf      { path "07.Metagene/" }
    // tss_minus_pdf     { path "07.Metagene/" }
    // tes_plus_pdf      { path "07.Metagene/" }
    // tes_minus_pdf     { path "07.Metagene/" }
}
