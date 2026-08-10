#!/usr/bin/env nextflow
//
// SUBWORKFLOW: deliver
// Chains: BigWig -> Pol II profiling -> IGV tracks -> GEO prep -> Report
//

include { bigwig_gen      } from '../modules/bigwig_gen.nf'
include { pol2_profiling  } from '../modules/pol2_profiling.nf'
include { igv_tracks      } from '../modules/igv_tracks.nf'
include { geo_prep        } from '../modules/geo_prep.nf'
include { report          } from '../modules/report.nf'

workflow deliver {
    take:
    strand_bams       // channel: tuple val(meta), path(plus_bam), path(minus_bam)
    bam_ch            // channel: tuple val(meta), path(bam)
    gene_body_counts  // channel: path — merged count matrix
    sample_sheet_ch   // channel: val(path)
    config_ch         // channel: val(config)

    main:
    // BigWig generation (per-sample, strand-specific)
    bigwig_gen(strand_bams)

    // Pol II active site profiling (per-sample)
    pol2_profiling(bam_ch, Channel.empty())

    // Collect BAMs and BigWigs for IGV
    all_bams = bam_ch.map { meta, bam -> bam }.collect()
    all_bws  = bigwig_gen.out.bigwigs
        .map { meta, pbw, mbw -> [pbw, mbw] }
        .flatten()
        .collect()

    // IGV track organization
    igv_tracks(all_bams, all_bws)

    // GEO submission prep
    geo_prep(sample_sheet_ch, gene_body_counts, config_ch)

    // Final report
    report(Channel.value(file("${launchDir}/results")))

    emit:
    bigwigs         = bigwig_gen.out.bigwigs
    profiling_xlsx  = pol2_profiling.out.profiling_xlsx
    tss_heatmap     = pol2_profiling.out.tss_heatmap
    igv_output      = igv_tracks.out.igv_output
    geo_output      = geo_prep.out.geo_output
    report_pdf      = report.out.report_pdf
}
