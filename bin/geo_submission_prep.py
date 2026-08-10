#!/usr/bin/env python3
"""
geo_submission_prep.py - Prepare GEO submission template files for PRO-seq data.

Generates:
    1. metadata.xlsx       - Sample metadata spreadsheet
    2. processed_data/      - Directory with count matrices, normalized bigWig symlinks
    3. raw_data_md5.txt     - MD5 checksums for raw data files

Usage:
    python geo_submission_prep.py --samplesheet samples.csv --counts gene_counts.txt \\
        --output_dir geo_submission [--bigwig_dir <dir>] [--raw_dir <dir>]
"""

import argparse
import csv
import hashlib
import logging
import os
import shutil
import sys
from datetime import datetime

try:
    import openpyxl
    from openpyxl.styles import Font, PatternFill, Alignment, Border, Side
    from openpyxl.utils import get_column_letter
    HAS_OPENPYXL = True
except ImportError:
    HAS_OPENPYXL = False
    print("WARNING: openpyxl not installed. Install with: pip install openpyxl")
    print("         Excel metadata file will be written as CSV instead.")

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
# GEO metadata template columns
# ---------------------------------------------------------------------------
GEO_METADATA_FIELDS = [
    "Sample title",
    "Organism",
    "Molecule type",
    "Library strategy",
    "Library source",
    "Library selection",
    "Library layout",
    "Platform",
    "Instrument model",
    "Description",
    "Data processing step",
    "Processed data file",
    "Raw file",
]

GEO_SERIES_FIELDS = [
    "Series title",
    "Series summary",
    "Series overall design",
    "Contributor",
    "Submission date",
    "Release date",
]


# ---------------------------------------------------------------------------
# Read samplesheet CSV
# Expected columns (minimum): sample, condition, replicate, organism, molecule, layout, platform, instrument
# Optional columns: description, raw_file, any additional metadata
# ---------------------------------------------------------------------------
def read_samplesheet(samplesheet_path):
    """Read the samplesheet and return a list of dictionaries."""
    if not os.path.exists(samplesheet_path):
        log.error("Samplesheet not found: %s", samplesheet_path)
        sys.exit(1)

    samples = []
    with open(samplesheet_path, "r", newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            # Strip whitespace from all values
            cleaned = {k.strip(): v.strip() if v else "" for k, v in row.items()}
            samples.append(cleaned)

    if not samples:
        log.error("No samples found in samplesheet")
        sys.exit(1)

    log.info("Read %d samples from %s", len(samples), samplesheet_path)
    log.info("  Columns: %s", ", ".join(samples[0].keys()))
    return samples


# ---------------------------------------------------------------------------
# Generate metadata rows from samples
# ---------------------------------------------------------------------------
def build_metadata_rows(samples, counts_file):
    """Build GEO-compliant metadata rows."""
    rows = []

    for samp in samples:
        sample_name = samp.get("sample", samp.get("sample_name", "unknown"))
        condition = samp.get("condition", samp.get("group", "unknown"))
        replicate = samp.get("replicate", "1")

        sample_title = f"PRO-seq {sample_name} ({condition}, rep{replicate})"

        organism = samp.get("organism", "Homo sapiens")
        molecule = samp.get("molecule", "polyA-depleted RNA")
        strategy = "PRO-seq"
        source = "TRANSCRIPTOMIC"
        selection = "cDNA"
        layout = samp.get("layout", "SINGLE")
        platform = samp.get("platform", "ILLUMINA")
        instrument = samp.get("instrument", samp.get("instrument_model", "NextSeq 2000"))
        description = samp.get("description", f"PRO-seq library: {condition}, replicate {replicate}")
        raw_file = samp.get("raw_file", f"{sample_name}.fastq.gz")
        processed_file = samp.get("processed_data", f"{sample_name}_plus.bw, {sample_name}_minus.bw")

        data_processing = (
            "Base calling: bcl2fastq v2.20; "
            "Adapter trimming: cutadapt v4.0; "
            "Alignment: STAR v2.7.11b; "
            "Filtering: samtools view -q 30; "
            "Deduplication: picard MarkDuplicates; "
            "BigWig generation: deepTools bamCoverage --normalizeUsing RPKM; "
            "Count matrix: featureCounts"
        )

        row = {
            "Sample title": sample_title,
            "Organism": organism,
            "Molecule type": molecule,
            "Library strategy": strategy,
            "Library source": source,
            "Library selection": selection,
            "Library layout": layout,
            "Platform": platform,
            "Instrument model": instrument,
            "Description": description,
            "Data processing step": data_processing,
            "Processed data file": processed_file,
            "Raw file": raw_file,
        }
        rows.append(row)

    log.info("Built %d metadata rows", len(rows))
    return rows


# ---------------------------------------------------------------------------
# Write metadata to Excel (or CSV fallback)
# ---------------------------------------------------------------------------
def write_metadata_xlsx(rows, output_path):
    """Write GEO metadata to Excel file."""
    if not rows:
        log.error("No metadata rows to write")
        return

    if HAS_OPENPYXL:
        _write_xlsx(rows, output_path)
    else:
        csv_path = output_path.replace(".xlsx", ".csv")
        _write_csv_fallback(rows, csv_path)
        log.warning("openpyxl not available; metadata written as CSV: %s", csv_path)


def _write_xlsx(rows, output_path):
    """Write metadata to .xlsx with formatting."""
    wb = openpyxl.Workbook()

    # --- Sheet 1: Metadata Template ---
    ws_meta = wb.active
    ws_meta.title = "METADATA_TEMPLATE"

    # Header style
    header_font = Font(name="Arial", size=11, bold=True, color="FFFFFF")
    header_fill = PatternFill(start_color="4472C4", end_color="4472C4", fill_type="solid")
    header_align = Alignment(horizontal="center", vertical="center", wrap_text=True)
    thin_border = Border(
        left=Side(style="thin"),
        right=Side(style="thin"),
        top=Side(style="thin"),
        bottom=Side(style="thin"),
    )

    # GEO instruction row
    ws_meta.merge_cells("A1:L1")
    ws_meta["A1"] = "GEO Metadata Template — PRO-seq Experiment"
    ws_meta["A1"].font = Font(name="Arial", size=14, bold=True)
    ws_meta.merge_cells("A2:L2")
    ws_meta["A2"] = (
        f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')} | "
        f"Samples: {len(rows)} | "
        "Fill in all fields and upload to GEO at https://www.ncbi.nlm.nih.gov/geo/info/seq.html"
    )
    ws_meta["A2"].font = Font(name="Arial", size=10, italic=True)

    # Write headers
    fields = GEO_METADATA_FIELDS
    for col_idx, field in enumerate(fields, 1):
        cell = ws_meta.cell(row=4, column=col_idx, value=field)
        cell.font = header_font
        cell.fill = header_fill
        cell.alignment = header_align
        cell.border = thin_border

    # Write data rows
    for row_idx, row_data in enumerate(rows, 5):
        for col_idx, field in enumerate(fields, 1):
            cell = ws_meta.cell(row=row_idx, column=col_idx, value=row_data.get(field, ""))
            cell.alignment = Alignment(wrap_text=True, vertical="top")
            cell.border = thin_border

    # Auto-fit column widths (approximate)
    for col_idx in range(1, len(fields) + 1):
        max_len = len(str(fields[col_idx - 1]))
        for row_idx in range(5, 5 + len(rows)):
            val = ws_meta.cell(row=row_idx, column=col_idx).value or ""
            max_len = max(max_len, min(len(str(val)), 60))
        ws_meta.column_dimensions[get_column_letter(col_idx)].width = max_len + 4

    # Freeze header
    ws_meta.freeze_panes = "A5"

    # --- Sheet 2: README ---
    ws_readme = wb.create_sheet("README")
    instructions = [
        "GEO Submission Instructions for PRO-seq Data",
        "",
        "1. Review and complete the METADATA_TEMPLATE sheet.",
        "2. Ensure all required fields are filled for each sample:",
        "   - Sample title",
        "   - Organism",
        "   - Molecule type (e.g., polyA-depleted RNA)",
        "   - Library strategy (PRO-seq)",
        "   - Library source (TRANSCRIPTOMIC)",
        "   - Library selection (cDNA)",
        "   - Library layout (SINGLE or PAIRED)",
        "   - Platform (ILLUMINA)",
        "   - Instrument model",
        "   - Description",
        "   - Data processing step",
        "3. Prepare processed data files in the processed_data/ directory.",
        "4. Prepare raw FASTQ files and note their MD5 checksums in raw_data_md5.txt.",
        "5. Upload all files to the GEO submission portal:",
        "   https://www.ncbi.nlm.nih.gov/geo/info/seq.html",
        "",
        "Common GEO field values for PRO-seq:",
        "  Molecule type: polyA-depleted RNA",
        "  Library strategy: PRO-seq (select 'OTHER' in GEO, specify 'PRO-seq')",
        "  Library source: TRANSCRIPTOMIC",
        "  Library selection: cDNA",
        "",
        "For paired-end PRO-seq (4sU + No-4sU):",
        "  Submit each sample individually.",
        "  Include Input (No-4sU) and Nascent (4sU-enriched) as separate samples.",
        "  Note the relationship in the sample descriptions.",
    ]
    for i, line in enumerate(instructions, 1):
        cell = ws_readme.cell(row=i, column=1, value=line)
        if i == 1:
            cell.font = Font(name="Arial", size=12, bold=True)
        elif line.startswith("  "):
            cell.font = Font(name="Consolas", size=10)

    ws_readme.column_dimensions["A"].width = 100

    wb.save(output_path)
    log.info("Metadata Excel written: %s", output_path)


def _write_csv_fallback(rows, output_path):
    """Fallback: write metadata as CSV."""
    fields = GEO_METADATA_FIELDS
    with open(output_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    log.info("Metadata CSV written: %s", output_path)


# ---------------------------------------------------------------------------
# Copy processed data files (count matrices, etc.)
# ---------------------------------------------------------------------------
def copy_processed_data(counts_file, bigwig_dir, output_dir):
    """Copy or symlink processed data files into the GEO submission directory."""
    processed_dir = os.path.join(output_dir, "processed_data")
    os.makedirs(processed_dir, exist_ok=True)

    copied = []

    # Copy count matrix
    if counts_file and os.path.exists(counts_file):
        dest = os.path.join(processed_dir, os.path.basename(counts_file))
        shutil.copy2(counts_file, dest)
        copied.append(dest)
        log.info("Copied count matrix: %s -> %s", counts_file, dest)
    else:
        log.warning("Counts file not found: %s", counts_file)

    # Copy/symlink bigWig files
    if bigwig_dir and os.path.isdir(bigwig_dir):
        for fname in os.listdir(bigwig_dir):
            if fname.endswith(".bw") or fname.endswith(".bigwig") or fname.endswith(".bigWig"):
                src = os.path.join(bigwig_dir, fname)
                dest = os.path.join(processed_dir, fname)
                if not os.path.exists(dest):
                    shutil.copy2(src, dest)
                    copied.append(dest)
        log.info("Copied %d bigWig files from %s", len(copied), bigwig_dir)
    elif bigwig_dir:
        log.warning("BigWig directory not found: %s", bigwig_dir)

    return copied


# ---------------------------------------------------------------------------
# Compute MD5 checksums
# ---------------------------------------------------------------------------
def compute_md5_checksums(raw_dir, output_path):
    """Compute MD5 checksums for all files in the raw data directory."""
    if not raw_dir or not os.path.isdir(raw_dir):
        log.warning("Raw data directory not found; skipping MD5 checksums: %s", raw_dir)
        return

    md5_entries = []
    for root, dirs, files in os.walk(raw_dir):
        for fname in sorted(files):
            fpath = os.path.join(root, fname)
            log.info("Computing MD5: %s", fname)
            md5_hash = _md5_file(fpath)
            rel_path = os.path.relpath(fpath, raw_dir)
            md5_entries.append(f"{md5_hash}  {rel_path}")

    with open(output_path, "w") as f:
        f.write("# MD5 checksums for raw PRO-seq data files\n")
        f.write(f"# Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"# Source directory: {raw_dir}\n")
        f.write("\n")
        f.write("\n".join(md5_entries))
        f.write("\n")

    log.info("MD5 checksums written (%d files): %s", len(md5_entries), output_path)


def _md5_file(filepath, chunk_size=8192):
    """Compute MD5 hash of a file."""
    md5 = hashlib.md5()
    with open(filepath, "rb") as f:
        while True:
            chunk = f.read(chunk_size)
            if not chunk:
                break
            md5.update(chunk)
    return md5.hexdigest()


# ---------------------------------------------------------------------------
# Generate series-level metadata summary
# ---------------------------------------------------------------------------
def write_series_metadata(output_dir, samples):
    """Write minimal series-level metadata placeholder."""
    series_path = os.path.join(output_dir, "series_metadata.txt")
    with open(series_path, "w") as f:
        f.write("GEO Series Metadata (fill in before submission)\n")
        f.write("=" * 60 + "\n\n")
        f.write("Series title: PRO-seq profiling of ...\n")
        f.write("Series summary: ...\n")
        f.write("Series overall design: ...\n")
        f.write(f"Contributor: ...\n")
        f.write(f"Submission date: {datetime.now().strftime('%Y-%m-%d')}\n")
        f.write("Release date: [specify release date]\n")
        f.write("\n")
        f.write("Sample summary:\n")
        for s in samples:
            f.write(f"  - {s.get('sample', 'unknown')}: {s.get('condition', '?')} "
                    f"(rep{s.get('replicate', '?')})\n")

    log.info("Series metadata placeholder: %s", series_path)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        description="Prepare GEO submission template files for PRO-seq data"
    )
    parser.add_argument(
        "--samplesheet",
        required=True,
        help="CSV samplesheet with columns: sample, condition, replicate, organism, layout, platform, instrument",
    )
    parser.add_argument(
        "--counts",
        default=None,
        help="Path to gene-level count matrix (tab-separated)",
    )
    parser.add_argument(
        "--output_dir",
        required=True,
        help="Output directory for GEO submission files",
    )
    parser.add_argument(
        "--bigwig_dir",
        default=None,
        help="Directory containing bigWig (.bw) files to copy into processed_data/",
    )
    parser.add_argument(
        "--raw_dir",
        default=None,
        help="Directory containing raw FASTQ files for MD5 checksum computation",
    )

    args = parser.parse_args()

    log.info("=== GEO Submission Preparation ===")

    os.makedirs(args.output_dir, exist_ok=True)

    # Read samplesheet
    samples = read_samplesheet(args.samplesheet)

    # Build metadata
    metadata_rows = build_metadata_rows(samples, args.counts)

    # Write metadata.xlsx
    metadata_path = os.path.join(args.output_dir, "metadata.xlsx")
    write_metadata_xlsx(metadata_rows, metadata_path)

    # Copy processed data
    copy_processed_data(args.counts, args.bigwig_dir, args.output_dir)

    # Compute MD5 checksums for raw data
    if args.raw_dir:
        md5_path = os.path.join(args.output_dir, "raw_data_md5.txt")
        compute_md5_checksums(args.raw_dir, md5_path)

    # Series-level metadata placeholder
    write_series_metadata(args.output_dir, samples)

    log.info("=== GEO Submission Package ===")
    log.info("Output directory: %s", os.path.abspath(args.output_dir))
    log.info("Files:")
    for item in sorted(os.listdir(args.output_dir)):
        item_path = os.path.join(args.output_dir, item)
        if os.path.isfile(item_path):
            log.info("  - %s", item)
        elif os.path.isdir(item_path):
            n_files = sum(1 for _ in os.listdir(item_path) if os.path.isfile(os.path.join(item_path, _)))
            log.info("  - %s/ (%d files)", item, n_files)

    log.info("Done. Review and complete templates before GEO submission.")


if __name__ == "__main__":
    main()
