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
//   4. Pausing index (TSS / gene body) + boxplot
//   5. Differential expression (DESeq2 on gene body counts)
//   6. Pol II active-site single-base distribution (bedtools) + bigWigs (deepTools)
//   7. Metagene TSS/TES profiles (deepTools)
//
// Samplesheet: sample,group,r1,r2
// ============================================================

include { parse_config      } from './modules/parse_config.nf'
include { preprocess        } from './subworkflows/preprocess.nf'
include { align_bowtie2     } from './subworkflows/align_bowtie2.nf'

include { prepare_genome    } from './subworkflows/prepare_genome.nf'
include { quantification    } from './subworkflows/quantification.nf'
include { pause_analysis    } from './subworkflows/pause_analysis.nf'
include { tss_meta          } from './subworkflows/tss_meta.nf'
include { diff              } from './subworkflows/diff.nf'
include { POL2_FIVEPRIME    } from './modules/bedtools/pol2_fiveprime.nf'

workflow {

    main:
    // ========================================================================
    // Step 0: Parse genome configuration
    // ========================================================================
    config_ch = parse_config().map { cfg_file ->
        def config = [:]
        cfg_file.text.split('\n').each { line ->
            def p = line.split(': ', 2)
            if (p.size() >= 2) { config[p[0].trim()] = p[1].trim() }
        }
        println "============================================"
        println " PRO-seq pipeline — config parsed"
        println " build           : ${config.build}"
        println " genome_fasta    : ${config.genome_fasta}"
        println " bowtie2_index   : ${config.bowtie2_index}"
        println " gtf             : ${config.gtf}"
        println " gene_annotation : ${config.gene_annotation}"
        println " sample_sheet    : ${params.sample_sheet}"
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
            [meta, [row.r1, row.r2]]
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
    // Step 4: quantification (per-sample featureCounts -> merged matrices)
    //   - one featureCounts run per sample per region (handles mixed SE/PE)
    //   - quantification merges per-sample counts into tss/genebody matrices
    // ========================================================================
    quantification(align_bowtie2.out.bam, prepare_genome.out.tss_saf, prepare_genome.out.genebody_saf)
    

    // ========================================================================
    // Step 5: Pausing index (TSS / gene body) + boxplot
    //   Groups come from the samplesheet 'group' column.
    // ========================================================================
    // Re-read the samplesheet (a tiny local CSV) to derive group -> [samples].
    // Avoids multicasting read_ch with `into`, which the DSL2 compiler failed to
    // resolve as a channel operator ("Missing process or function into").
    groups_config_ch = Channel.fromPath(params.sample_sheet)
        .splitCsv(header: true)
        .map { row -> [row.group ?: 'unknown', row.sample] }
        .toList()
        .map { pairs ->
            def groups = [:]
            pairs.each { group, sample ->
                if (!groups.containsKey(group)) groups[group] = []
                groups[group] << sample
            }
            return groups
        }

    pause_analysis(quantification.out.tss_matrix, quantification.out.genebody_matrix, groups_config_ch)

    // ========================================================================
    // Step 5b: Differential expression (DESeq2 on gene body counts)
    //   groups come from the samplesheet; comparisons from params.compared_groups;
    //   annotation from the database config (gene_annotation).
    // ========================================================================
    annotation_ch = config_ch.map { it -> it.gene_annotation ?: '' }

    diff(quantification.out.genebody_matrix, groups_config_ch, annotation_ch)

    // ========================================================================
    // Step 5c: Pol II active-site single-base distribution (5' end coverage)
    //   Per sample: signed + / - bedGraphs (bedtools genomecov).
    //   bigWig coverage is now produced inside TSS_meta (Step 7).
    // ========================================================================
    POL2_FIVEPRIME(align_bowtie2.out.bam)

    // ========================================================================
    // Step 7: Metagene (TSS profile via deepTools, per-sample)
    // ========================================================================
    gtf_ch = config_ch.map { it -> it.gtf }
    tss_meta(align_bowtie2.out.bam, align_bowtie2.out.bai, gtf_ch)

    // ========================================================================
    // Publish results to output directories
    // ========================================================================
    publish:
    fastqc_zip         = preprocess.out.fastqc_zip
    fastqc_html        = preprocess.out.fastqc_html
    fastp_json         = preprocess.out.fastp_json
    fastp_html         = preprocess.out.fastp_html
    fastp_log          = preprocess.out.fastp_log
    trimmed_reads      = preprocess.out.trimmed_reads
    base_quality_plot  = preprocess.out.base_quality_plot
    statistics         = preprocess.out.statistics
    statistics_log     = preprocess.out.statistics_log

    bam                = align_bowtie2.out.bam
    bai                = align_bowtie2.out.bai
    alignRate          = align_bowtie2.out.alignRate

    tss_counts         = quantification.out.tss_counts
    genebody_counts    = quantification.out.genebody_counts
    tss_matrix         = quantification.out.tss_matrix
    genebody_matrix    = quantification.out.genebody_matrix

    pi_all             = pause_analysis.out.pi_all
    pi_boxplot         = pause_analysis.out.pi_boxplot

    diff_results       = diff.out.results

    pol2_plus_bedgraph   = POL2_FIVEPRIME.out.plus_bedgraph
    pol2_minus_bedgraph  = POL2_FIVEPRIME.out.minus_bedgraph
    pol2_signed_bedgraph = POL2_FIVEPRIME.out.signed_bedgraph

    coverage_bw         = tss_meta.out.bigwig

    tss_meta_matrix     = tss_meta.out.matrix
    tss_meta_profile    = tss_meta.out.profile
}

// ------------------------------------------------------------------
// Output directive
// ------------------------------------------------------------------
output {
    fastqc_zip        { path "02.fastqc/" }
    fastqc_html       { path "02.fastqc/" }
    fastp_json        { path "03.Data_QC/" }
    fastp_html        { path "03.Data_QC/" }
    fastp_log         { path "03.Data_QC/" }
    trimmed_reads     { path "03.Data_QC/" }
    base_quality_plot { path "03.Data_QC/" }
    statistics        { path "03.Data_QC/" }
    statistics_log    { path "03.Data_QC/" }

    bam               { path "04.Alignment/" }
    bai               { path "04.Alignment/" }
    alignRate         { path "04.Alignment/" }

    tss_counts        { path "05.Quantification/" }
    genebody_counts   { path "05.Quantification/" }
    tss_matrix        { path "05.Quantification/" }
    genebody_matrix   { path "05.Quantification/" }

    pi_all            { path "06.Pausing_Index/" }
    pi_boxplot        { path "06.Pausing_Index/" }

    diff_results      { path "07.Differential_Expression/" }

    pol2_plus_bedgraph   { path "08.Pol2_Active_Site/" }
    pol2_minus_bedgraph  { path "08.Pol2_Active_Site/" }
    pol2_signed_bedgraph { path "08.Pol2_Active_Site/" }

    coverage_bw         { path "09.DeepTools_Coverage/" }

    tss_meta_matrix     { path "10.TSS_Metagene/" }
    tss_meta_profile    { path "10.TSS_Metagene/" }
}
