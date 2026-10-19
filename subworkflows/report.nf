#!/usr/bin/env nextflow
//
// SUBWORKFLOW: report
// 报告打包链（params.report 独立开关，与 compared_groups 无关）：
//   1. OVERALL_STAT:  read_statistics + bowtie2 比对率 -> overall_statistics.txt
//   2. TXT2XLSX:      txt2xlsx.R ×5 -> report_xlsx/（5 张交付 xlsx）
//   3. PAUSING_XLSX:  Pausing_Index.tsv + promoter.bed/genebody.bed + 注释表 -> PROSeq_pausing.xlsx
//   4. REPORT:        Report.R 打包 -> 16 章交付目录
// 无差异分析时 diff/enrich/plot 通道为空：TXT2XLSX 收空 diff_dir（全 [skip]），
// REPORT 收空 List 建空 stage 目录，Report.R 对应章节只 [skip] 不失败。
//

include { OVERALL_STAT     } from '../modules/report/overall_stat.nf'
include { TXT2XLSX         } from '../modules/report/txt2xlsx.nf'
include { PAUSING_XLSX     } from '../modules/report/pausing_xlsx.nf'
include { POL2_SIGNAL_XLSX } from '../modules/report/pol2_signal_xlsx.nf'
include { REPORT           } from '../modules/report/report.nf'

workflow report {
    take:
    read_statistics      // path: read_statistics.txt（单条）
    align_rate           // tuple(meta, *_summary_bowtie2.txt) × 样本
    base_quality_plot    // path: *.tiff/*.pdf × 样本（过滤前/后碱基质量图）
    diff_dir             // path: deseq2_out/（0..1 条；未做差异分析时为空通道）
    profile_file         // path: genebody.annotated.txt（normalize.R 产物，main_mapped 口径表达谱）
    enrich_files         // path glob：gokegg/gsea 子目录（0..n）
    de_plot_files        // path glob：heatmap/scatter/volcano（0..n）
    metagene_files       // path：metagene plot + matrix（n 条，必非空）
    group_heatmaps       // path：组级 metagene heatmap（metagene_group.out.plot；可空）
    pol2_signal          // path：pol2_signal_table.tsv（signal_table 产物；full 模式为空通道）
    pol2_signal_note     // path：pol2_signal_table.note.txt（signal_table 产物；full 模式为空通道）
    pi_table             // path：Pausing_Index.tsv（单条）
    promoter_bed         // path：promoter.bed（单条，代表 transcript TSS 窗口）
    genebody_bed         // path：genebody.bed（单条，代表 transcript gene body）
    annotation           // val：gene_annotation 路径字符串（可为 ''）
    config_yml           // path: params.yml
    samplesheet_csv      // path: samplesheet.csv
    resolved_config      // path: config.txt

    main:
    bowtie2_summaries = align_rate.map { _meta, f -> f }.collect()
    OVERALL_STAT(read_statistics, bowtie2_summaries)

    // 无差异分析时 diff_dir 为空通道 → 传空串（进程内 [ -n ] 全 skip）
    diff_dir_path = diff_dir.map { d -> d.toString() }.ifEmpty('')
    TXT2XLSX(read_statistics, OVERALL_STAT.out.overall_statistics, diff_dir_path, profile_file)

    // signal_table note → txt2xlsx 转 PROSeq_pol2_signal.xlsx（full 模式空通道 → 进程跳过、xlsx 空）
    POL2_SIGNAL_XLSX(pol2_signal, pol2_signal_note)

    PAUSING_XLSX(pi_table, promoter_bed, genebody_bed, annotation)

    REPORT(
        enrich_files.collect().ifEmpty([]),      // 空 → 单个空 List 元素，val 绑定合法
        metagene_files.collect(),
        PAUSING_XLSX.out.xlsx.collect(),
        de_plot_files.collect().ifEmpty([]),
        base_quality_plot.collect(),
        group_heatmaps.collect().ifEmpty([]),
        POL2_SIGNAL_XLSX.out.xlsx.map { f -> f.toString() }.ifEmpty(''),
        TXT2XLSX.out.xlsx_dir,
        config_yml,
        samplesheet_csv,
        resolved_config
    )

    emit:
    report = REPORT.out.report
}
