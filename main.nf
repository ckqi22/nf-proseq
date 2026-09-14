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
//   9. Report packaging（Report.R 打包交付；params.report 开关，独立于差异分析）
//
// ============================================================

include { parse_config                } from './modules/parse_config.nf'
include { prepare_genome              } from './subworkflows/prepare_genome.nf'
include { preprocess                  } from './subworkflows/preprocess.nf'
include { align_bowtie2               } from './subworkflows/align_bowtie2.nf'
include { spikein                     } from './subworkflows/spikein.nf'
include { quantification              } from './subworkflows/quantification.nf'
include { diff                        } from './subworkflows/diff.nf'
include { enrich                      } from './subworkflows/enrich.nf'
include { pol2_count                  } from './subworkflows/pol2_count.nf'
include { pause_analysis              } from './subworkflows/pause_analysis.nf'
include { metagene                    } from './subworkflows/metagene.nf'
include { metagene_group              } from './subworkflows/metagene_group.nf'
include { profile as profile_genebody } from './subworkflows/profile.nf'
include { profile as profile_promoter } from './subworkflows/profile.nf'
include { SIGNAL_TABLE                } from './modules/signal_table.nf'
include { report as report_package    } from './subworkflows/report.nf'

workflow {

    main:
    // ========================================================================
    // Step 0: Parse genome configuration -> Read samplesheet
    // ========================================================================
    parse_config_raw = parse_config()
    
    config_ch = parse_config_raw.map { cfg_file ->
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
        println " spike_fasta     : ${config.spike_fasta}"
        println " spike_index     : ${config.spike_index}"
        println " spike_chroms    : ${config.spike_chroms}"
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
    // Spike-in（条件分支）：params.spike_genome 配置时，统计每样本 spike reads，
    // 把 spike_count 并入 meta 供 GENOMECOV 缩放 bigWig，并生成因子表供 profile 表缩放；
    // 未配置时照旧（库大小 CPM），spike_factors_ch 给空串。
    // ========================================================================
    spike_enabled = params.spike_genome?.trim() || params.spike_fasta?.trim() || params.spike_index?.trim()

    if (spike_enabled) {
        spike = spikein(align_bowtie2.out.bam_bai, prepare_genome.out.spike_chroms)
        // spike_count 并入 meta（新 map，不改共享 meta，避免污染 quantification 消费的 meta）
        bam_for_pol2 = align_bowtie2.out.bam
            .join(spike.counts.map { m, f -> [m, f.text.trim().toInteger()] }, by: [0])
            .map { meta, bam, c -> [meta + [spike_count: c], bam] }
        spike_factors_ch = spike.factors.map { it.toRealPath().toString() }
        spike_factors_out = spike.factors
    } else {
        bam_for_pol2 = align_bowtie2.out.bam
        spike_factors_ch = channel.value('')
        spike_factors_out = channel.empty()
    }

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
        diff_result_ch    = diff.out.result
        checkde_result_ch = diff.out.checkde_result
        de_plot_ch        = diff.out.de_plot
        pca_plot_ch       = diff.out.pca_plot
    } else {
        diff_result_ch    = channel.empty()
        checkde_result_ch = channel.empty()
        de_plot_ch        = channel.empty()
        pca_plot_ch       = channel.empty()
    }

    // ========================================================================
    // Step 5: Enrichment analysis
    // ========================================================================
    if (params.compared_groups) {
        enrich(diff.out.passed)
        enrich_result_ch = enrich.out.gokegg_result.mix(enrich.out.gsea_result)
    } else {
        enrich_result_ch = channel.empty()
    }

    // ========================================================================
    // Step 6: Pol II active-site single-base distribution
    //   (produces per-base bedGraph/bigWig + single-base promoter/genebody matrices for PI)
    // ========================================================================
    pol2_count(bam_for_pol2, prepare_genome.out.promoter_bed, prepare_genome.out.genebody_bed)

    // Step 6b: 生成带注释 + 原始 count + 标准化值的 profile 表
    methods = params.normalize_methods ?: 'cpm,fpkm'
    profile_genebody(quantification.out.genebody_matrix, annotation_ch, methods, 'genebody', spike_factors_ch)
    // profile_promoter(pol2_count.out.promoter_matrix, annotation_ch, methods, 'promoter')

    // ========================================================================
    // Step 7: TSS Metagene profiles（deepTools：每样本 profile + heatmap；
    //   组图按 samplesheet 的 group 列，仅 ≥2 样本的显式分组；'unknown'/单样本组不出组图）
    // ========================================================================
    // signal_mode：single=单碱基5'端(默认) / full=full read 全长覆盖度（metagene 用 CPM 版）
    signal_mode = (params.signal_mode?.trim() ?: 'single')

    // 同一 channel 直接喂两个子流程（DSL2 多消费者，无需 into 分叉）
    metagene_bw = (signal_mode == 'full') ? pol2_count.out.bigwig_coverage_cpm
                                          : pol2_count.out.bigwig_cpm
    metagene(metagene_bw, prepare_genome.out.tss_bed)
    metagene_group(metagene_bw, prepare_genome.out.tss_bed)

    // 报告用混合通道：publish 块原先在此处 mix 会独占两个源通道，
    // 上移到这里定义一次，publish 直接引用、报告侧再 mix（各通道保持单一算子消费者）
    tss_metagene_plot_all   = metagene.out.plot.mix(metagene_group.out.plot)
    tss_metagene_matrix_all = metagene.out.matrix.mix(metagene_group.out.matrix)

    // ========================================================================
    // Step 8: Pause index (single-base promoter / gene body)
    // ========================================================================
    pause_analysis(pol2_count.out.promoter_matrix, pol2_count.out.genebody_matrix, groups_config_ch)

    // ========================================================================
    // Step 9: 报告打包（params.report 独立开关；与 compared_groups 无关，
    //         无差异分析也打包，diff/enrich/plot 相关章节由 Report.R [skip]）
    // ========================================================================
    if (params.report) {
        report_package(
            preprocess.out.statistics,
            align_bowtie2.out.alignRate,
            preprocess.out.fastqc_raw_zip,
            preprocess.out.fastqc_trimmed_zip,
            preprocess.out.base_quality_plot,
            diff_result_ch,
            enrich_result_ch,
            de_plot_ch,
            tss_metagene_plot_all.mix(tss_metagene_matrix_all),
            metagene_group.out.plot,
            pause_analysis.out.pi_all,
            prepare_genome.out.promoter_bed.map { _name, file -> file },
            prepare_genome.out.genebody_bed.map { _name, file -> file },
            annotation_ch,
            file(params.report_config),
            file(params.sample_sheet),
            parse_config_raw
        )
        report_ch = report_package.out.report
    } else {
        report_ch = channel.empty()
    }

    // ========================================================================
    // Step 10（暂不接入）: 逐碱基 Pol II 活性位点信号表（按组聚合 + RPM + gene/transcript 注释）
    // ========================================================================
    // bg_files    = pol2_count.out.bedGraph
    //                 .map { _meta, plus, minus -> [plus, minus] }
    //                 .collect()
    // bg_manifest = pol2_count.out.bedGraph
    //                 .map { meta, plus, minus -> "${meta.sample}\t${plus.name}\t${minus.name}" }
    //                 .collect()
    //                 .map { lines -> lines.join('\n') }

    // groups_yml_ch = groups_config_ch.map { groups ->
    //     def lines = []
    //     groups.each { name, samples ->
    //         lines << "${name}:"
    //         samples.each { s -> lines << "  - ${s}" }
    //     }
    //     lines.join('\n')
    // }

    // SIGNAL_TABLE(bg_files, bg_manifest, groups_yml_ch,
    //              prepare_genome.out.gene_bed, config_ch.map { it.gtf }, annotation_ch)

    // ========================================================================
    // Publish results to output directories
    // ========================================================================
    publish:
    info                = parse_config.out.info
    representative_gtf  = prepare_genome.out.representative_gtf
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

    promoter_pol2_counts    = pol2_count.out.promoter_counts.map { _meta, file -> file }
    genebody_pol2_counts    = pol2_count.out.genebody_counts.map { _meta, file -> file }
    genebody_counts         = quantification.out.genebody_counts
    promoter_pol2_matrix    = pol2_count.out.promoter_matrix
    genebody_pol2_matrix    = pol2_count.out.genebody_matrix
    genebody_matrix         = quantification.out.genebody_matrix

    diff_result         = diff_result_ch
    checkde_result      = checkde_result_ch
    de_plot             = de_plot_ch
    pca_plot            = pca_plot_ch

    enrich_result       = enrich_result_ch

    coverage_bw         = pol2_count.out.bigwig
    coverage_bw_cpm     = pol2_count.out.bigwig_cpm
    coverage_full_bw    = pol2_count.out.bigwig_coverage_cpm
    // pol2_signal_table   = SIGNAL_TABLE.out.signal_table

    genebody_profile    = profile_genebody.out.annotated
    // promoter_profile    = profile_promoter.out.annotated

    spike_factors       = spike_factors_out

    tss_metagene_plot   = tss_metagene_plot_all
    tss_metagene_matrix = tss_metagene_matrix_all
    group_avg_bigwig    = metagene_group.out.avg_bigwig
    metagene_combined_plot   = metagene.out.combined_plot.mix(metagene_group.out.combined_plot)
    metagene_combined_matrix = metagene.out.combined_matrix.mix(metagene_group.out.combined_matrix)

    pi_all              = pause_analysis.out.pi_all
    pi_boxplot          = pause_analysis.out.pi_boxplot
    report              = report_ch

}

// ------------------------------------------------------------------
// Output directive
// ------------------------------------------------------------------
output {
    info                { path "01.Info/" }
    representative_gtf  { path "01.Info/" }
    promoter_bed        { path "01.Info/" }
    genebody_bed        { path "01.Info/" }
    genebody_union_saf  { path "01.Info/" }
    gene_bed            { path "01.Info/" }

    fastqc_raw_zip      { path "02.Fastqc/" }
    fastqc_raw_html     { path "02.Fastqc/" }
    fastqc_trimmed_zip  { path "02.Fastqc/" }
    fastqc_trimmed_html { path "02.Fastqc/" }

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
    genebody_matrix         { path "05.Quantification/" }
    genebody_profile        { path "05.Quantification/" }
    // promoter_profile        { path "05.Quantification/" }
    spike_factors           { path "05.Quantification/" }

    diff_result         { path "06.Differential_Expression/" }
    checkde_result      { path "06.Differential_Expression/" }
    de_plot             { path "06.Differential_Expression/" }
    pca_plot            { path "06.Differential_Expression/" }

    enrich_result       { path "07.enrich/" }

    coverage_bw         { path "08.Pol2_coverage/" }
    coverage_bw_cpm     { path "08.Pol2_coverage/" }
    coverage_full_bw    { path "08.Pol2_coverage/" }
    // pol2_signal_table   { path "08.Pol2_coverage/" }

    tss_metagene_plot   {path "09.TSS_Metagene/"}
    tss_metagene_matrix {path "09.TSS_Metagene/"}
    group_avg_bigwig    { path "09.TSS_Metagene/" }
    metagene_combined_plot   { path "09.TSS_Metagene/" }
    metagene_combined_matrix { path "09.TSS_Metagene/" }

    pi_all              { path "10.Pausing_Index/" }
    pi_boxplot          { path "10.Pausing_Index/" }

    report              { path "11.Report/" }
}
