#!/usr/bin/env python3
"""
Download Cell Painting Gallery images for CCM2 analysis.

Downloads TIFF images for specific wells from AWS S3 Cell Painting Gallery.

Usage:
    pixi run download-images --wells A19,I21,M19,C23 --sites 1,2,3

    # Download all CCM2 wells (ref: A19, I21; var: M19, C23)
    pixi run download-images
"""
import argparse
import subprocess
import os
from pathlib import Path
from typing import List

# AWS S3 bucket
S3_BUCKET = "s3://cellpainting-gallery/cpg0052-jump-varchamp"
S3_SOURCE = "source_1"

# CCM2 well positions
CCM2_WELLS = {
    "ref": ["A19", "I21"],      # Reference (CCM2)
    "var": ["M19", "C23"],      # Variant (CCM2_Val53Ile)
}

# Batch plate mappings
BATCHES = {
    "2025_06_10_Batch_18": {
        "2025_06_02_B18A8A10R1_P2T1": "2025_06_02_B18A8A10R1_P2T1__2025-06-02T08_19_26-Measurement1",
        "2025_06_02_B18A8A10R1_P2T2": "2025_06_02_B18A8A10R1_P2T2__2025-06-02T09_29_48-Measurement1",
        "2025_06_02_B18A8A10R1_P2T3": "2025_06_02_B18A8A10R1_P2T3__2025-06-02T11_08_44-Measurement1",
        "2025_06_02_B18A8A10R1_P2T4": "2025_06_02_B18A8A10R1_P2T4__2025-06-02T12_18_34-Measurement1",
    },
    "2025_06_10_Batch_19": {
        "2025_06_03_B19A8A10R1_P2T1": "2025_06_03_B19A8A10R1_P2T1__2025-06-03T13_27_13-Measurement1",
        "2025_06_04_B19A8A10R1_P2T2": "2025_06_04_B19A8A10R1_P2T2__2025-06-04T08_17_36-Measurement2",
        "2025_06_04_B19A8A10R1_P2T3": "2025_06_04_B19A8A10R1_P2T3__2025-06-04T09_30_37-Measurement1",
        "2025_06_04_B19A8A10R1_P2T4": "2025_06_04_B19A8A10R1_P2T4__2025-06-04T11_49_06-Measurement1",
    },
}

# Channel mapping
CHANNELS = {
    "DAPI": 1,
    "GFP": 2,
    "AGP": 3,
    "Mito": 4,
}

# Row letter to number
LETTER_TO_ROW = {chr(i): i - 64 for i in range(65, 81)}

# Output directory
OUTPUT_DIR = Path(__file__).parent.parent / "data/raw/images"


def well_to_rowcol(well: str) -> tuple:
    """Convert well position to (row, col) numbers."""
    row = LETTER_TO_ROW[well[0].upper()]
    col = int(well[1:])
    return row, col


def build_image_filename(row: int, col: int, site: int, channel: int) -> str:
    """Build TIFF filename."""
    return f"r{row:02d}c{col:02d}f{site:02d}p01-ch{channel}sk1fk1fl1.tiff"


def download_well_images(
    batch: str,
    plate: str,
    plate_dir: str,
    well: str,
    sites: List[int],
    channels: List[str],
    output_dir: Path,
    dry_run: bool = False
) -> int:
    """
    Download images for a single well.

    Returns number of files downloaded.
    """
    row, col = well_to_rowcol(well)

    # Build S3 path
    s3_base = f"{S3_BUCKET}/{S3_SOURCE}/{batch}/images/{plate_dir}/Images"

    # Local output path
    local_base = output_dir / batch / plate / well
    local_base.mkdir(parents=True, exist_ok=True)

    downloaded = 0
    for site in sites:
        for ch_name in channels:
            ch_num = CHANNELS[ch_name]
            filename = build_image_filename(row, col, site, ch_num)

            s3_path = f"{s3_base}/{filename}"
            local_path = local_base / filename

            if local_path.exists():
                print(f"  Exists: {local_path.name}")
                continue

            if dry_run:
                print(f"  Would download: {filename}")
                continue

            # Download using aws s3 cp
            cmd = ["aws", "s3", "cp", "--no-sign-request", s3_path, str(local_path)]
            try:
                subprocess.run(cmd, check=True, capture_output=True)
                print(f"  Downloaded: {filename}")
                downloaded += 1
            except subprocess.CalledProcessError as e:
                print(f"  Failed: {filename} - {e.stderr.decode()[:100]}")

    return downloaded


def main():
    parser = argparse.ArgumentParser(description="Download CPG images for CCM2 analysis")
    parser.add_argument(
        "--wells",
        default=",".join(CCM2_WELLS["ref"] + CCM2_WELLS["var"]),
        help="Wells to download (comma-separated, default: A19,I21,M19,C23)"
    )
    parser.add_argument(
        "--sites",
        default="1,2,3,4,5",
        help="Sites to download (comma-separated, default: 1,2,3,4,5)"
    )
    parser.add_argument(
        "--channels",
        default="DAPI,GFP,AGP,Mito",
        help="Channels to download (default: all four)"
    )
    parser.add_argument(
        "--batches",
        default="18,19",
        help="Batches to download (18, 19, or both)"
    )
    parser.add_argument(
        "--plates",
        default=None,
        help="Specific plates (default: all plates in batches)"
    )
    parser.add_argument(
        "--output-dir",
        default=str(OUTPUT_DIR),
        help="Output directory for images"
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Show what would be downloaded without downloading"
    )

    args = parser.parse_args()

    wells = [w.strip() for w in args.wells.split(",")]
    sites = [int(s.strip()) for s in args.sites.split(",")]
    channels = [c.strip() for c in args.channels.split(",")]
    batch_nums = [b.strip() for b in args.batches.split(",")]
    output_dir = Path(args.output_dir)

    print("=" * 60)
    print("CCM2 Cell Painting Image Download")
    print("=" * 60)
    print(f"Wells: {wells}")
    print(f"Sites: {sites}")
    print(f"Channels: {channels}")
    print(f"Output: {output_dir}")
    if args.dry_run:
        print("DRY RUN - no files will be downloaded")
    print()

    total_downloaded = 0

    for batch_num in batch_nums:
        batch_id = f"2025_06_10_Batch_{batch_num}"
        if batch_id not in BATCHES:
            print(f"Unknown batch: {batch_num}")
            continue

        print(f"\nBatch {batch_num}:")
        print("-" * 40)

        plates = BATCHES[batch_id]
        if args.plates:
            plate_filter = [p.strip() for p in args.plates.split(",")]
            plates = {k: v for k, v in plates.items() if k in plate_filter}

        for plate, plate_dir in plates.items():
            print(f"\n  Plate: {plate}")

            for well in wells:
                print(f"    Well {well}:")
                n = download_well_images(
                    batch_id, plate, plate_dir, well,
                    sites, channels, output_dir, args.dry_run
                )
                total_downloaded += n

    print()
    print("=" * 60)
    if args.dry_run:
        print("DRY RUN complete")
    else:
        print(f"Downloaded {total_downloaded} files")
    print("=" * 60)


if __name__ == "__main__":
    main()
