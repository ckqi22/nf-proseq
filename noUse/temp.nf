include { strand_split  } from '../modules/strand_split.nf'
    strand_split(alignment_bowtie2.out.bam)

    // ------------------------------------------------------------------
    // Step 4: Strand separation
    // Forward reads (-F 0x10) go to plus.bam for plus-strand genes
    // Reverse reads (-f 0x10) go to minus.bam for minus-strand genes
    // ------------------------------------------------------------------

    // Strand indices
    plus_bai  = strand_split.out.plus_bai
    minus_bai = strand_split.out.minus_bai

    emit:
    // Strand-split BAMs: tuple val(meta), path(plus_bam), path(minus_bam)
    strand_bams = strand_split.out.strand_bams
