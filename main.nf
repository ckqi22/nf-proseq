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
include { REMOVE_RRNA                 } from './modules/bowtie2/remove_rrna.nf'
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
    // Step 1b: 去 rRNA（opt-in，比对前）
    //   bowtie2 比对 rRNA_index，只按 R1 判定，保留 R1 未比对 rRNA 的 reads；
    //   remove_rrna=false 或物种段未配 rRNA_index 时直接用 trimmed reads。
    // ========================================================================
    rrna_index_ch = config_ch.map { it -> it.rRNA_index ?: '' }
                             .filter { s -> s.trim() }
                             .map { idx -> [[id: 'rrna'], file(idx).parent] }

    if (params.remove_rrna) {
        rrna_removed    = REMOVE_RRNA(preprocess.out.trimmed_reads, rrna_index_ch)
        reads_for_align = rrna_removed.clean_reads
        rrna_rate_ch    = rrna_removed.rrna_rate
    } else {
        reads_for_align = preprocess.out.trimmed_reads
        rrna_rate_ch    = channel.empty()
    }

    // ========================================================================
    // Step 2: alignment
    // ========================================================================
    spike_enabled = params.spike_genome?.trim() || params.spike_fasta?.trim() || params.spike_index?.trim()

    // spike_chroms 名单透传给 EXTRACT_R1：开 spike 时按名单剔 spike 染色体出纯主 r1_bam；
    // 关 spike 时空串占位（EXTRACT_R1 不剔，SE 走 cp 快路径）。
    spike_chroms_for_align = spike_enabled
        ? prepare_genome.out.spike_chroms.map { it -> it.toString() }
        : channel.value('')

    align_bowtie2(reads_for_align, prepare_genome.out.index, prepare_genome.out.fasta, spike_chroms_for_align)

    // ========================================================================
    // Spike-in（条件分支）：spike_enabled 时统计 spike reads，spike_count 并入 meta 供缩放；
    // 未配置时按库大小 CPM。total_mapped（read1 mapped，含 spike）恒并入 meta，供 spike 占比质控；
    // main_mapped（read1 mapped，纯主，EXTRACT_R1 剔 spike 后）供算 scale_cpm。
    // ========================================================================

    // 计数文件 → 整数通道（total_mapped / main_mapped 恒有；spike_count 仅开 spike 时有）
    def to_int = { m, f -> [m, f.text.trim().toInteger()] }

    total_mapped_ch = align_bowtie2.out.total_mapped.map(to_int)
    main_mapped_ch  = align_bowtie2.out.main_mapped.map(to_int)

    if (spike_enabled) {
        spike = spikein(align_bowtie2.out.bam_bai, prepare_genome.out.spike_chroms, align_bowtie2.out.total_mapped)
        spike_count_ch = spike.counts.map(to_int)
        spike_factors_ch = spike.factors                          // path(spikein_scale_factors.tsv)
        spike_factors_out = spike.factors
    } else {
        spike_count_ch = channel.empty()
        spike_factors_ch = channel.value([])                      // 空 list → NORMALIZE 判空（不加 --spike_factors）
        spike_factors_out = channel.empty()
    }

    // join 都在「原始 meta」上做（join by:[0] 比较整个 meta，先写任何额外 key 会失配）。
    // spike 关时 spike_count_ch 空，remainder:true 左外连 → sc=null，不并入 spike_count。
    // main_mapped 由 EXTRACT_R1 直接产出（纯主，已剔 spike），不再 tm - sc 反推。
    bam_for_pol2 = align_bowtie2.out.r1_bam
        .join(total_mapped_ch, by: [0])
        .join(spike_count_ch, by: [0], remainder: true)
        .join(main_mapped_ch, by: [0])
        .map { meta, r1_bam, tm, sc, mm ->
            def extra = [total_mapped: tm, main_mapped: mm]
            if (sc != null) extra.spike_count = sc
            [meta + extra, r1_bam]
        }

    // ========================================================================
    // Step 3: quantification (per-sample featureCounts -> merged gene body matrix)
    // ========================================================================
    quantification(align_bowtie2.out.r1_bam, prepare_genome.out.genebody_union_saf)
    
    // ========================================================================
    // Step 4: Differential expression (DESeq2 on gene body counts)
    // ========================================================================
    annotation_ch = config_ch.map { it -> it.gene_annotation ?: '' }

    if (params.diff?.compared_groups) {
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
    if (params.diff?.compared_groups) {
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
    profile_genebody(quantification.out.genebody_matrix, annotation_ch, methods, 'genebody', spike_factors_ch,
                     align_bowtie2.out.main_mapped.map { _m, f -> f }.collect())
    // profile_promoter(pol2_count.out.promoter_matrix, annotation_ch, methods, 'promoter')

    // ========================================================================
    // Step 7: TSS Metagene profiles（deepTools：每样本 profile + heatmap；
    //   组图按 samplesheet 的 group 列，仅 ≥2 样本的显式分组；'unknown'/单样本组不出组图）
    // ========================================================================
    // signal_mode：single=单碱基活性位点端(reverse 库=R1 5'端 / forward 库=R1 3'端，见 params.strandedness)
    //              full=full read 全长覆盖度
    //              both=两者都算(对比用)
    signal_mode = (params.signal_mode?.trim() ?: 'single')

    // 下游 metagene/metagene_group 的分析 bigWig：归一化优先 spike（开 spike 用 _spike，否则 _cpm），
    // signal_type 取 full/single（both 落单碱基）。spike_enabled 是 main 作用域布尔（上方已算，非 channel）。
    metagene_bw = (signal_mode == 'full')
        ? (spike_enabled ? pol2_count.out.bigwig_full_spike : pol2_count.out.bigwig_full_cpm)
        : (spike_enabled ? pol2_count.out.bigwig_spike      : pol2_count.out.bigwig_cpm)
    metagene(metagene_bw, prepare_genome.out.tss_bed)
    metagene_group(metagene_bw, prepare_genome.out.tss_bed, prepare_genome.out.chrom_sizes)

    // 报告用混合通道：publish 块原先在此处 mix 会独占两个源通道，
    // 上移到这里定义一次，publish 直接引用、报告侧再 mix（各通道保持单一算子消费者）
    tss_metagene_plot_all   = metagene.out.plot.mix(metagene_group.out.plot)
    tss_metagene_matrix_all = metagene.out.matrix.mix(metagene_group.out.matrix)

    // ========================================================================
    // Step 8: Pause index (single-base promoter / gene body)
    // ========================================================================
    pause_analysis(pol2_count.out.promoter_matrix, pol2_count.out.genebody_matrix, groups_config_ch)


    // ========================================================================
    // Step 9: 逐碱基 Pol II 活性位点信号表（活跃基因 pause 窗口，raw + 归一化双轨）
    //   manifest 6 列：group \t sample \t plus_raw \t minus_raw \t plus_norm \t minus_norm（bigWig basename）
    //   归一化轨：开 spike 用 spike 版，否则 CPM（与 metagene_bw 同口径）；raw 轨恒为单碱基原始计数。
    //   需单碱基活性位点端信号（signal_mode=single/both；末端由 params.strandedness 决定：reverse→5'端 / forward→3'端）；
    //   full 模式无单碱基 bigWig，跳过。
    // ========================================================================
    if (signal_mode != 'full') {
        norm_bigwig   = spike_enabled ? pol2_count.out.bigwig_spike : pol2_count.out.bigwig_cpm
        raw_norm_bigwig = pol2_count.out.bigwig.join(norm_bigwig, by: [0])

        manifest = raw_norm_bigwig
            .map { meta, plus_raw, minus_raw, plus_norm, minus_norm -> [meta.group, meta.sample, plus_raw.name, minus_raw.name, plus_norm.name, minus_norm.name].join('\t') }
            .collect()
            .map { lines -> lines.join('\n') }

        raw_norm_bws = raw_norm_bigwig
            .map { _meta, plus_raw, minus_raw, plus_norm, minus_norm -> [plus_raw, minus_raw, plus_norm, minus_norm] }
            .flatten()
            .collect()

        SIGNAL_TABLE(
            manifest,
            raw_norm_bws,
            prepare_genome.out.promoter_bed.map { _name, f -> f },
            pol2_count.out.promoter_matrix,
            pol2_count.out.genebody_matrix,
            prepare_genome.out.representative_gtf,
            annotation_ch,
            params.strandedness?.trim() ?: 'reverse'
        )
        signal_table_ch = SIGNAL_TABLE.out.signal_table
        signal_table_note_ch = SIGNAL_TABLE.out.note
    } else {
        signal_table_ch = channel.empty()      // full 模式无单碱基 bigWig，跳过
        signal_table_note_ch = channel.empty()
    }


    // ========================================================================
    // Step 10: 报告打包（params.report 独立开关；与 compared_groups 无关，
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
            signal_table_ch,
            signal_table_note_ch,
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
    // Publish results to output directories
    // ========================================================================
    publish:
    info                = parse_config.out.info
    representative_gtf  = prepare_genome.out.representative_gtf
    promoter_bed        = prepare_genome.out.promoter_bed.map { _name, file -> file }
    genebody_bed        = prepare_genome.out.genebody_bed.map { _name, file -> file }
    genebody_union_saf  = prepare_genome.out.genebody_union_saf.map { _name, file -> file }
    gene_bed            = prepare_genome.out.gene_bed
    tss_bed             = prepare_genome.out.tss_bed
    chrom_sizes         = prepare_genome.out.chrom_sizes
    spike_chroms        = prepare_genome.out.spike_chroms

    fastqc_raw_zip      = preprocess.out.fastqc_raw_zip.map { _meta, file -> file }
    fastqc_raw_html     = preprocess.out.fastqc_raw_html.map { _meta, file -> file }
    fastqc_trimmed_zip  = preprocess.out.fastqc_trimmed_zip.map { _meta, file -> file }
    fastqc_trimmed_html = preprocess.out.fastqc_trimmed_html.map { _meta, file -> file }
    fastp_json          = preprocess.out.fastp_json.map { _meta, file -> file }
    fastp_html          = preprocess.out.fastp_html
    fastp_log           = preprocess.out.fastp_log
    trimmed_reads       = preprocess.out.trimmed_reads.map { _meta, file -> file }
    base_quality_plot   = preprocess.out.base_quality_plot
    statistics          = preprocess.out.statistics

    bam                 = align_bowtie2.out.bam.map { _meta, file -> file }
    bai                 = align_bowtie2.out.bai.map { _meta, file -> file }
    unmapped_bam        = align_bowtie2.out.unmapped_bam.map { _meta, file -> file }
    alignRate           = align_bowtie2.out.alignRate.map { _meta, file -> file }
    rrna_rate           = rrna_rate_ch.map { _meta, file -> file }

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

    coverage_bw             = pol2_count.out.bigwig.map { _meta, plus, minus -> [plus, minus] }
    coverage_bw_cpm         = pol2_count.out.bigwig_cpm.map { _meta, plus, minus -> [plus, minus] }
    coverage_bw_spike       = pol2_count.out.bigwig_spike.map { _meta, plus, minus -> [plus, minus] }
    coverage_full_bw        = pol2_count.out.bigwig_full.map { _meta, plus, minus -> [plus, minus] }
    coverage_full_bw_cpm    = pol2_count.out.bigwig_full_cpm.map { _meta, plus, minus -> [plus, minus] }
    coverage_full_bw_spike  = pol2_count.out.bigwig_full_spike.map { _meta, plus, minus -> [plus, minus] }
    pol2_signal_table       = signal_table_ch
    genebody_profile        = profile_genebody.out.annotated
    // promoter_profile    = profile_promoter.out.annotated

    spike_factors               = spike_factors_out
    metagene_plot               = tss_metagene_plot_all
    metagene_matrix             = tss_metagene_matrix_all
    group_bw                    = metagene_group.out.avg_bigwig
    metagene_combined_plot      = metagene.out.combined_plot.mix(metagene_group.out.combined_plot)
    metagene_combined_matrix    = metagene.out.combined_matrix.mix(metagene_group.out.combined_matrix)

    pi_all      = pause_analysis.out.pi_all
    pi_boxplot  = pause_analysis.out.pi_boxplot

    report      = report_ch
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
    tss_bed             { path "01.Info/" }
    chrom_sizes         { path "01.Info/" }
    spike_chroms        { path "01.Info/" }

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
    unmapped_bam        { path "04.Alignment/" }
    alignRate           { path "04.Alignment/" }
    rrna_rate           { path "04.Alignment/" }

    promoter_pol2_counts    { path "05.Quantification/" }
    genebody_pol2_counts    { path "05.Quantification/" }
    genebody_counts         { path "05.Quantification/" }
    promoter_pol2_matrix    { path "05.Quantification/" }
    genebody_pol2_matrix    { path "05.Quantification/" }
    genebody_matrix         { path "05.Quantification/" }
    genebody_profile        { path "05.Quantification/" }
    spike_factors           { path "05.Quantification/" }
    // promoter_profile        { path "05.Quantification/" }

    diff_result     { path "06.Differential_Expression/" }
    checkde_result  { path "06.Differential_Expression/" }
    de_plot         { path "06.Differential_Expression/" }
    pca_plot        { path "06.Differential_Expression/" }

    enrich_result   { path "07.enrich/" }

    coverage_bw             { path "08.Pol2_coverage/" }
    coverage_bw_cpm         { path "08.Pol2_coverage/" }
    coverage_bw_spike       { path "08.Pol2_coverage/" }
    coverage_full_bw        { path "08.Pol2_coverage/" }
    coverage_full_bw_cpm    { path "08.Pol2_coverage/" }
    coverage_full_bw_spike  { path "08.Pol2_coverage/" }
    group_bw                { path "08.Pol2_coverage/" }
    pol2_signal_table       { path "08.Pol2_coverage/" }

    metagene_plot               { path "09.TSS_Metagene/" }
    metagene_matrix             { path "09.TSS_Metagene/" }
    metagene_combined_plot      { path "09.TSS_Metagene/" }
    metagene_combined_matrix    { path "09.TSS_Metagene/" }

    pi_all      { path "10.Pausing_Index/" }
    pi_boxplot  { path "10.Pausing_Index/" }

    report  { path "11.Report/" }
}
