#!/usr/bin/env python3
"""
Extract CCM2 and CCM2_Val53Ile profiles from VarChAMP Batch 18-19.

Creates two types of profiles:
1. Feature-selected GFP profiles (for standard classification)
2. Full GFP profiles with IntegratedIntensity (for GFP-corrected classification)

Usage:
    pixi run extract-data
"""
import polars as pl
from pathlib import Path

# Source paths
VARCHAMP_DIR = Path("/data/users/shenrunx/igvf/varchamp/2025_varchamp_snakemake")
BATCH_PROFILES_DIR = VARCHAMP_DIR / "2.snakemake_pipeline/outputs/batch_profiles"

# Feature-selected profiles (for standard classification)
FEATSELECT_PIPELINE = "profiles_tcdropped_filtered_var_mad_outlier_featselect_filtcells.parquet"
# Pre-feature-selection profiles (for GFP-corrected with IntegratedIntensity)
FULL_PIPELINE = "profiles_tcdropped_filtered_var_mad_outlier.parquet"

# Output path
OUTPUT_DIR = Path(__file__).parent.parent / "data/raw/single_cell_profiles"

# Alleles to extract
ALLELES = ["CCM2", "CCM2_Val53Ile"]

# Batches
BATCHES = {
    "2025_06_10_Batch_18": "batch18",
    "2025_06_10_Batch_19": "batch19",
}


def get_gfp_columns(df: pl.LazyFrame) -> list:
    """Get GFP feature columns only."""
    cols = df.collect_schema().names()
    meta_cols = [c for c in cols if c.startswith("Metadata_")]
    gfp_cols = [c for c in cols if "GFP" in c and not c.startswith("Metadata_")]
    return meta_cols + gfp_cols


def create_cell_id(df: pl.DataFrame) -> pl.DataFrame:
    """Create unique CellID from plate, well, image, object."""
    return df.with_columns(
        pl.concat_str([
            pl.col("Metadata_Plate"),
            pl.col("Metadata_Well"),
            pl.col("Metadata_ImageNumber").cast(pl.Utf8),
            pl.col("Metadata_ObjectNumber").cast(pl.Utf8),
        ], separator="_").alias("CellID")
    )


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    for batch_id, batch_name in BATCHES.items():
        print(f"\n{'='*60}")
        print(f"Processing {batch_id}")
        print("=" * 60)

        # =====================================================================
        # 1. Feature-selected GFP profiles (for standard classification)
        # =====================================================================
        featselect_path = BATCH_PROFILES_DIR / batch_id / FEATSELECT_PIPELINE
        out_gfp_path = OUTPUT_DIR / f"{batch_name}_ccm2_profiles.parquet"

        print(f"\n1. Loading feature-selected profiles: {featselect_path.name}")
        df_featselect = pl.scan_parquet(featselect_path)

        # Get GFP columns only
        gfp_cols = get_gfp_columns(df_featselect)
        print(f"   Selecting {len(gfp_cols)} columns (GFP + metadata)")

        # Filter to CCM2 alleles and select GFP columns
        df_gfp = (
            df_featselect
            .select(gfp_cols)
            .filter(pl.col("Metadata_gene_allele").is_in(ALLELES))
            .collect()
        )
        df_gfp = create_cell_id(df_gfp)

        print(f"   CCM2 cells: {df_gfp.shape[0]}")
        print(f"   Alleles: {df_gfp['Metadata_gene_allele'].unique().to_list()}")

        # Save
        df_gfp.write_parquet(out_gfp_path)
        size_mb = out_gfp_path.stat().st_size / 1024 / 1024
        print(f"   Saved to: {out_gfp_path.name} ({size_mb:.1f} MB)")

        # =====================================================================
        # 2. Full GFP profiles with IntegratedIntensity (for GFP-corrected)
        # =====================================================================
        full_path = BATCH_PROFILES_DIR / batch_id / FULL_PIPELINE
        out_full_gfp_path = OUTPUT_DIR / f"{batch_name}_ccm2_profiles_full_gfp.parquet"

        print(f"\n2. Loading full profiles: {full_path.name}")
        df_full = pl.scan_parquet(full_path)

        # Get GFP columns only
        full_gfp_cols = get_gfp_columns(df_full)
        print(f"   Selecting {len(full_gfp_cols)} columns (GFP + metadata)")

        # Create CellID for filtering
        df_full_with_id = (
            df_full
            .select(full_gfp_cols)
            .with_columns(
                pl.concat_str([
                    pl.col("Metadata_Plate"),
                    pl.col("Metadata_Well"),
                    pl.col("Metadata_ImageNumber").cast(pl.Utf8),
                    pl.col("Metadata_ObjectNumber").cast(pl.Utf8),
                ], separator="_").alias("CellID")
            )
        )

        # Filter to only cells that passed QC (in feature-selected profiles)
        valid_cell_ids = df_gfp["CellID"].to_list()
        print(f"   Filtering to {len(valid_cell_ids)} valid cells...")

        df_full_filtered = (
            df_full_with_id
            .filter(pl.col("CellID").is_in(valid_cell_ids))
            .collect()
        )

        print(f"   Matched cells: {df_full_filtered.shape[0]}")

        # Verify IntegratedIntensity column exists
        int_cols = [c for c in df_full_filtered.columns if "IntegratedIntensity" in c and "GFP" in c]
        print(f"   IntegratedIntensity columns: {int_cols[:3]}...")

        # Save
        df_full_filtered.write_parquet(out_full_gfp_path)
        size_mb = out_full_gfp_path.stat().st_size / 1024 / 1024
        print(f"   Saved to: {out_full_gfp_path.name} ({size_mb:.1f} MB)")

    print(f"\n{'='*60}")
    print("Done!")
    print("=" * 60)
    print("\nOutput files:")
    for f in sorted(OUTPUT_DIR.glob("*.parquet")):
        size_mb = f.stat().st_size / 1024 / 1024
        print(f"  {f.name}: {size_mb:.1f} MB")


if __name__ == "__main__":
    main()
