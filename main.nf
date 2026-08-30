#!/usr/bin/env nextflow
//
// ============================================================
// nf_PROseq — Nextflow DSL2 PRO-seq Analysis Pipeline
// ============================================================
//
// Workflow:
//   1. QC: Fastp QC + adapter trimming
//   2. Alignment: Bowtie2 / ...
//   3. Quantification: promoter / genebody
//   4. Differential expression (DESeq2 on gene body counts)
//   5. Enrichment analysis (GO/KEGG/GSEA)
//   6. Pol II active-site single-base distribution
//   7. Metagene TSS profiles
//   8. Pausing index
//
// ============================================================

include { parse_config      } from './modules/parse_config.nf'
include { prepare_genome    } from './subworkflows/prepare_genome.nf'
include { preprocess        } from './subworkflows/preprocess.nf'
include { align_bowtie2     } from './subworkflows/align_bowtie2.nf'
include { quantification    } from './subworkflows/quantification.nf'
include { diff              } from './subworkflows/diff.nf'
include { enrich            } from './subworkflows/enrich.nf'
include { pol2_count      } from './subworkflows/pol2_count.nf'
include { pause_analysis    } from './subworkflows/pause_analysis.nf'
include { tss_meta          } from './subworkflows/tss_meta.nf'
include { profile as profile_genebody } from './subworkflows/profile.nf'
include { profile as profile_promoter } from './subworkflows/profile.nf'
include { SIGNAL_TABLE        } from './modules/signal_table.nf'

workflow {

    main:
    // ========================================================================
    // Step 0: Parse genome configuration -> Read samplesheet
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
        println "============================================"
        return config
    }

    // Resolve reference genome (fasta + bowtie2 index) and prepare
    // promoter + genebody SAF annotations from the reference GTF
    prepare_genome(config_ch)

    // ========================================================================
    // Read samplesheet
    // ========================================================================
    read_ch = channel.fromPath(params.sample_sheet)
        .splitCsv(header: true)
        .map { row ->
            if (!row.sample) { error "samplesheet missing 'sample' column" }
            if (!row.r1)     { error "samplesheet missing 'r1' column" }
            def meta = [
                sample:     row.sample,
                group:      row.group ?: 'unknown',
                single_end: row.r2 ? false : true
            ]
            def reads = meta.single_end ? [row.r1] : [row.r1, row.r2]
            [meta, reads]
        }

    groups_config_ch = channel.fromPath(params.sample_sheet)
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

    // ========================================================================
    // Step 1: preprocess
    // ========================================================================
    preprocess(read_ch, params.adapter ?: 'I')

    // ========================================================================
    // Step 2: alignment
    // ========================================================================
    align_bowtie2(preprocess.out.trimmed_reads, prepare_genome.out.index, prepare_genome.out.fasta)
        
    // ========================================================================
    // Step 3: quantification (per-sample featureCounts -> merged gene body matrix)
    // ========================================================================
    quantification(align_bowtie2.out.bam, prepare_genome.out.genebody_union_saf)
    
    // ========================================================================
    // Step 4: Differential expression (DESeq2 on gene body counts)
    // ========================================================================
    annotation_ch = config_ch.map { it -> it.gene_annotation ?: '' }

    if (params.compared_groups) {
        diff(quantification.out.genebody_matrix, groups_config_ch, annotation_ch)
        diff_result_ch = diff.out.result
    } else {
        diff_result_ch = channel.empty()
    }

    // ========================================================================
    // Step 5: Enrichment analysis
    // ========================================================================
    if (params.compared_groups) {
        enrich(diff_result_ch)
        enrich_result_ch = enrich.out.gokegg_result.mix(enrich.out.gsea_result)
    } else {
        enrich_result_ch = channel.empty()
    }

    // ========================================================================
    // Step 6: Pol II active-site single-base distribution
    //   (produces per-base bedGraph/bigWig + single-base promoter/genebody matrices for PI)
    // ========================================================================
    pol2_count(align_bowtie2.out.bam, prepare_genome.out.promoter_bed, prepare_genome.out.genebody_bed)

    // Step 6b: 生成带注释 + 原始 count + 标准化值的 profile 表
    methods = params.normalize_methods ?: 'cpm,fpkm'
    profile_genebody(quantification.out.genebody_matrix, annotation_ch, methods, 'genebody')
    // profile_promoter(pol2_count.out.promoter_matrix, annotation_ch, methods, 'promoter')

    // ========================================================================
    // Step 7: Metagene TSS profiles
    // ========================================================================
    // gtf_ch = config_ch.map { it -> it.gtf }
    tss_meta(pol2_count.out.bigwig, prepare_genome.out.gene_bed)

    // ========================================================================
    // Step 8: Pause index (single-base promoter / gene body)
    // ========================================================================
    pause_analysis(pol2_count.out.promoter_matrix, pol2_count.out.genebody_matrix, groups_config_ch)

    // ========================================================================
    // Step 9: 逐碱基 Pol II 活性位点信号表（按组聚合 + RPM + gene/transcript 注释）
    // ========================================================================
    pol2_count.out.bedgraph
        .collect()
        .multiMap { tuples ->
            files:    tuples.collectMany { _meta, plus, minus -> [plus, minus] }
            manifest: tuples.collect { meta, plus, minus -> "${meta.sample}\t${plus.name}\t${minus.name}" }.join('\n')
        }
        .set { bg_mm }

    groups_yml_ch = groups_config_ch.map { groups ->
        def lines = []
        groups.each { name, samples ->
            lines << "${name}:"
            samples.each { s -> lines << "  - ${s}" }
        }
        lines.join('\n')
    }

    SIGNAL_TABLE(bg_mm.files, bg_mm.manifest, groups_yml_ch,
                 prepare_genome.out.gene_bed, config_ch.map { it.gtf }, annotation_ch)

    // ========================================================================
    // Publish results to output directories
    // ========================================================================
    publish:
    info                = parse_config.out.info
    promoter_bed        = prepare_genome.out.promoter_bed.map { _name, file -> file }
    genebody_bed        = prepare_genome.out.genebody_bed.map { _name, file -> file }
    genebody_union_saf  = prepare_genome.out.genebody_union_saf.map { _name, file -> file }
    gene_bed            = prepare_genome.out.gene_bed

    fastqc_raw_zip      = preprocess.out.fastqc_raw_zip
    fastqc_raw_html     = preprocess.out.fastqc_raw_html
    fastqc_trimmed_zip  = preprocess.out.fastqc_trimmed_zip
    fastqc_trimmed_html = preprocess.out.fastqc_trimmed_html
    fastp_json          = preprocess.out.fastp_json
    fastp_html          = preprocess.out.fastp_html
    fastp_log           = preprocess.out.fastp_log
    trimmed_reads       = preprocess.out.trimmed_reads
    base_quality_plot   = preprocess.out.base_quality_plot
    statistics          = preprocess.out.statistics

    bam                 = align_bowtie2.out.bam
    bai                 = align_bowtie2.out.bai
    alignRate           = align_bowtie2.out.alignRate

    promoter_pol2_counts    = pol2_count.out.promoter_counts
    genebody_pol2_counts    = pol2_count.out.genebody_counts
    genebody_counts         = quantification.out.genebody_counts
    promoter_pol2_matrix    = pol2_count.out.promoter_matrix
    genebody_pol2_matrix    = pol2_count.out.genebody_matrix
    genebody_matrix         = quantification.out.genebody_matrix

    diff_result         = diff_result_ch

    enrich_result       = enrich_result_ch

    coverage_bw         = pol2_count.out.bigwig
    pol2_signal_table   = SIGNAL_TABLE.out.signal_table

    genebody_profile    = profile_genebody.out.annotated
    // promoter_profile    = profile_promoter.out.annotated

    tss_meta_matrix     = tss_meta.out.matrix
    tss_meta_profile    = tss_meta.out.profile

    pi_all              = pause_analysis.out.pi_all
    pi_boxplot          = pause_analysis.out.pi_boxplot

}

// ------------------------------------------------------------------
// Output directive
// ------------------------------------------------------------------
output {
    info                { path "01.info/" }
    promoter_bed        { path "01.info/" }
    genebody_bed        { path "01.info/" }
    genebody_union_saf  { path "01.info/" }
    gene_bed            { path "01.info/" }

    fastqc_raw_zip      { path "02.fastqc/" }
    fastqc_raw_html     { path "02.fastqc/" }
    fastqc_trimmed_zip  { path "02.fastqc/" }
    fastqc_trimmed_html { path "02.fastqc/" }
    fastp_json          { path "03.Data_QC/" }
    fastp_html          { path "03.Data_QC/" }
    fastp_log           { path "03.Data_QC/" }
    trimmed_reads       { path "03.Data_QC/" }
    base_quality_plot   { path "03.Data_QC/" }
    statistics          { path "03.Data_QC/" }

    bam                 { path "04.Alignment/" }
    bai                 { path "04.Alignment/" }
    alignRate           { path "04.Alignment/" }

    promoter_pol2_counts    { path "05.Quantification/" }
    genebody_pol2_counts    { path "05.Quantification/" }
    genebody_counts         { path "05.Quantification/" }
    promoter_pol2_matrix    { path "05.Quantification/" }
    genebody_pol2_matrix    { path "05.Quantification/" }
    genebody_profile        { path "05.Quantification/" }
    // promoter_profile        { path "05.Quantification/" }

    diff_result         { path "06.Differential_Expression/" }

    enrich_result       { path "07.enrich/" }

    coverage_bw         { path "08.Pol2_coverage/" }
    pol2_signal_table   { path "08.Pol2_coverage/" }

    tss_meta_matrix     { path "09.TSS_Metagene/" }
    tss_meta_profile    { path "09.TSS_Metagene/" }

    pi_all              { path "10.Pausing_Index/" }
    pi_boxplot          { path "10.Pausing_Index/" }


}
