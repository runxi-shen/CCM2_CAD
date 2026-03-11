#!/usr/bin/env python3
"""
Extract cell location data for CCM2 visualization.

Creates a CSV file with cell metadata and pixel coordinates for cropping.
"""
import polars as pl
from pathlib import Path

# Source paths
VARCHAMP_DIR = Path("/data/users/shenrunx/igvf/varchamp/2025_varchamp_snakemake")
BATCH_PROFILES_DIR = VARCHAMP_DIR / "2.snakemake_pipeline/outputs/batch_profiles"
# Use pre-normalization profiles to get raw pixel coordinates
PROFILE_FILE = "profiles_tcdropped.parquet"

# Output path
OUTPUT_DIR = Path(__file__).parent.parent / "data/interim/cell_locations"

# Batches and CCM2 alleles
BATCHES = {
    "2025_06_10_Batch_18": "batch18",
    "2025_06_10_Batch_19": "batch19",
}
ALLELES = ["CCM2", "CCM2_Val53Ile"]

# Columns to extract
COLUMNS = [
    # Identifiers
    "Metadata_Plate",
    "Metadata_well_position",
    "Metadata_ImageNumber",
    "Metadata_ObjectNumber",
    "Metadata_gene_allele",
    "Metadata_node_type",
    # Cell locations (pixel coordinates)
    "Cells_Location_Center_X",
    "Cells_Location_Center_Y",
    "Nuclei_Location_Center_X",
    "Nuclei_Location_Center_Y",
]


def create_cell_id(df: pl.DataFrame) -> pl.DataFrame:
    """Create unique CellID from plate, well, image, object."""
    return df.with_columns(
        pl.concat_str([
            pl.col("Metadata_Plate"),
            pl.col("Metadata_well_position"),
            pl.col("Metadata_ImageNumber").cast(pl.Utf8),
            pl.col("Metadata_ObjectNumber").cast(pl.Utf8),
        ], separator="_").alias("CellID")
    )


def compute_site_number(df: pl.DataFrame) -> pl.DataFrame:
    """
    Compute site number (1-9) from ImageNumber.

    ImageNumber is a global ID. Each plate/well has 9 consecutive images (sites).
    Site = ImageNumber - min(ImageNumber for this plate/well) + 1
    """
    return df.with_columns(
        (pl.col("Metadata_ImageNumber") -
         pl.col("Metadata_ImageNumber").min().over(["Metadata_Plate", "Metadata_well_position"]) + 1
        ).alias("site")
    )


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    all_locations = []

    for batch_id, batch_name in BATCHES.items():
        print(f"\nProcessing {batch_id}...")

        profile_path = BATCH_PROFILES_DIR / batch_id / PROFILE_FILE
        if not profile_path.exists():
            print(f"  Profile not found: {profile_path}")
            continue

        # Load and filter
        df = (
            pl.scan_parquet(profile_path)
            .filter(pl.col("Metadata_gene_allele").is_in(ALLELES))
            .select(COLUMNS)
            .collect()
        )

        df = create_cell_id(df)
        df = compute_site_number(df)
        df = df.with_columns(pl.lit(batch_name).alias("batch"))

        print(f"  Cells: {len(df)}")
        print(f"  Alleles: {df['Metadata_gene_allele'].unique().to_list()}")

        all_locations.append(df)

    # Combine and save
    df_all = pl.concat(all_locations)

    # Save as CSV (small, easy to inspect)
    output_path = OUTPUT_DIR / "ccm2_cell_locations.csv"
    df_all.write_csv(output_path)

    # Also save as parquet for efficiency
    df_all.write_parquet(OUTPUT_DIR / "ccm2_cell_locations.parquet")

    print(f"\n{'='*60}")
    print(f"Total cells: {len(df_all)}")
    print(f"Saved to: {output_path}")
    print()
    print("Sample data:")
    print(df_all.head(5))


if __name__ == "__main__":
    main()
