#!/usr/bin/env python3
"""
bigwig_wrapper.py - Wrapper around deepTools bamCoverage for strand-specific bigwig generation.

For each BAM file in the input directory, generates:
  - <sample>_plus.bw  (forward strand)
  - <sample>_minus.bw (reverse strand)

Supports normalization methods: RPKM, CPM, BPM, or none.

Dependencies:
    deepTools (bamCoverage) must be installed and on PATH.
    samtools must be installed and on PATH (for BAM indexing).

Usage:
    python bigwig_wrapper.py --bam_dir <dir> --output_dir <dir> [--normalize RPKM] [--bin_size 10]
"""

import argparse
import os
import sys
import subprocess
import glob
import logging
from datetime import datetime

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------
logging.basicConfig(
    level=logging.INFO,
    format="[%(asctime)s] %(levelname)-5s %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Check required tools
# ---------------------------------------------------------------------------
def check_dependencies():
    """Ensure deepTools bamCoverage and samtools are available."""
    missing = []

    for tool in ["bamCoverage", "samtools"]:
        result = subprocess.run(
            ["which", tool], capture_output=True, text=True
        )
        if result.returncode != 0:
            missing.append(tool)

    if missing:
        log.error("Missing required tools: %s", ", ".join(missing))
        log.error("Install deepTools: conda install -c bioconda deeptools")
        log.error("Install samtools:  conda install -c bioconda samtools")
        sys.exit(1)

    log.info("deepTools bamCoverage and samtools found on PATH")


# ---------------------------------------------------------------------------
# Discover BAM files
# ---------------------------------------------------------------------------
def discover_bams(bam_dir):
    """Find all BAM files in the directory."""
    if not os.path.isdir(bam_dir):
        log.error("BAM directory not found: %s", bam_dir)
        sys.exit(1)

    bam_files = sorted(glob.glob(os.path.join(bam_dir, "*.bam")))
    bam_files += sorted(glob.glob(os.path.join(bam_dir, "*.BAM")))

    if not bam_files:
        log.error("No BAM files found in: %s", bam_dir)
        sys.exit(1)

    log.info("Found %d BAM files in %s", len(bam_files), bam_dir)
    for bf in bam_files:
        log.info("  %s", os.path.basename(bf))

    return bam_files


# ---------------------------------------------------------------------------
# Ensure BAM index exists
# ---------------------------------------------------------------------------
def ensure_index(bam_path):
    """Create BAM index (.bai) if missing."""
    bai_path = bam_path + ".bai"
    alt_bai = bam_path.replace(".bam", ".bai")

    if os.path.exists(bai_path) or os.path.exists(alt_bai):
        return

    log.info("Indexing: %s", os.path.basename(bam_path))
    result = subprocess.run(
        ["samtools", "index", bam_path],
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        log.error("samtools index failed for %s: %s", bam_path, result.stderr)
        sys.exit(1)


# ---------------------------------------------------------------------------
# Normalization argument mapping
# ---------------------------------------------------------------------------
NORMALIZE_MAP = {
    "RPKM": "RPKM",
    "CPM": "CPM",
    "BPM": "BPM",
    "NONE": None,
    "none": None,
}


# ---------------------------------------------------------------------------
# Generate bigwig for one BAM, one strand
# ---------------------------------------------------------------------------
def run_bamcoverage(
    bam_path, output_bw, strand, normalize, bin_size, effective_genome_size=None
):
    """
    Run deepTools bamCoverage for the specified strand.

    Parameters
    ----------
    bam_path : str
        Path to input BAM file.
    output_bw : str
        Path to output bigwig file.
    strand : str
        Strand filter: 'forward' for plus strand, 'reverse' for minus strand.
    normalize : str or None
        Normalization method (RPKM, CPM, BPM) or None for raw.
    bin_size : int
        Bin size in bp.
    effective_genome_size : int or None
        Effective genome size for RPKM normalization.
    """
    cmd = [
        "bamCoverage",
        "--bam", bam_path,
        "--outFileName", output_bw,
        "--outFileFormat", "bigwig",
        "--binSize", str(bin_size),
        "--filterRNAstrand", strand,
        "--numberOfProcessors", "max",
    ]

    if normalize:
        cmd.extend(["--normalizeUsing", normalize])
        if normalize == "RPKM" and effective_genome_size:
            cmd.extend(["--effectiveGenomeSize", str(effective_genome_size)])

    log.info("  Running: %s", " ".join(cmd))
    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        log.error("bamCoverage failed:")
        log.error("  stdout: %s", result.stdout.strip())
        log.error("  stderr: %s", result.stderr.strip())
        return False

    if os.path.exists(output_bw):
        size_kb = os.path.getsize(output_bw) / 1024
        log.info("  Created: %s (%.1f KB)", os.path.basename(output_bw), size_kb)
        return True
    else:
        log.error("  Output not created: %s", output_bw)
        return False


# ---------------------------------------------------------------------------
# Process a single BAM file
# ---------------------------------------------------------------------------
def process_bam(bam_path, output_dir, normalize, bin_size, effective_genome_size):
    """Generate plus and minus bigwig files for a single BAM."""
    sample_name = os.path.basename(bam_path).replace(".bam", "").replace(".BAM", "")

    ensure_index(bam_path)

    plus_bw  = os.path.join(output_dir, f"{sample_name}_plus.bw")
    minus_bw = os.path.join(output_dir, f"{sample_name}_minus.bw")

    log.info("Processing sample: %s", sample_name)

    success = True

    # Forward strand (plus)
    if not run_bamcoverage(
        bam_path, plus_bw, "forward", normalize, bin_size, effective_genome_size
    ):
        success = False

    # Reverse strand (minus)
    if not run_bamcoverage(
        bam_path, minus_bw, "reverse", normalize, bin_size, effective_genome_size
    ):
        success = False

    return success


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        description="Generate strand-specific bigWig files for PRO-seq BAM files using deepTools bamCoverage"
    )
    parser.add_argument(
        "--bam_dir",
        required=True,
        help="Directory containing input BAM files",
    )
    parser.add_argument(
        "--output_dir",
        required=True,
        help="Directory for output bigWig files",
    )
    parser.add_argument(
        "--normalize",
        default="RPKM",
        choices=["RPKM", "CPM", "BPM", "none"],
        help="Normalization method (default: RPKM). Use 'none' for raw counts.",
    )
    parser.add_argument(
        "--bin_size",
        type=int,
        default=10,
        help="Bin size in bp for bigWig summarization (default: 10)",
    )
    parser.add_argument(
        "--effective_genome_size",
        type=int,
        default=None,
        help="Effective genome size for RPKM normalization (e.g., 2913022398 for hg38)",
    )

    args = parser.parse_args()

    normalize = NORMALIZE_MAP.get(args.normalize, args.normalize)
    if args.normalize == "none":
        normalize = None

    log.info("=== PRO-seq BigWig Generation ===")
    log.info("BAM directory:     %s", args.bam_dir)
    log.info("Output directory:  %s", args.output_dir)
    log.info("Normalization:     %s", args.normalize if normalize else "None (raw counts)")
    log.info("Bin size:          %d bp", args.bin_size)

    # Setup
    check_dependencies()
    os.makedirs(args.output_dir, exist_ok=True)

    bam_files = discover_bams(args.bam_dir)

    # Process each BAM
    n_success = 0
    n_total = 0

    for bam_path in bam_files:
        n_total += 1
        if process_bam(bam_path, args.output_dir, normalize, args.bin_size, args.effective_genome_size):
            n_success += 1

    log.info("=== Summary ===")
    log.info("Total BAMs processed:  %d", n_total)
    log.info("Successful:            %d", n_success)
    log.info("Failed:                %d", n_total - n_success)

    if n_success < n_total:
        sys.exit(1)


if __name__ == "__main__":
    main()
