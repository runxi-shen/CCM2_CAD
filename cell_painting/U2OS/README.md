# CCM2_Val53Ile Cell Painting Analysis

XGBoost classification of CCM2_Val53Ile variant using GFP channel cell painting features from U2OS cells.

## Quick Start

```bash
# Install pixi if needed
curl -fsSL https://pixi.sh/install.sh | bash

# Set up environment (first time only)
cd /path/to/CCM2_CAD/cell_painting/U2OS
pixi install

# Run classification
pixi run classify
```

## Overview

This repository contains a standalone analysis of the CCM2_Val53Ile variant from VarChAMP Batches 18-19. The analysis:

## Data Source

Cell profiles are extracted from the VarChAMP pipeline: https://github.com/broadinstitute/2025_varchamp_snakemake

Two parquet file types are provided per batch:
- `batchNN_ccm2_profiles.parquet`: De-correlated, feature-selected GFP profiles (for standard classification)
- `batchNN_ccm2_profiles_full_gfp.parquet`: All normalized GFP features including IntegratedIntensity (for GFP-corrected analysis)

- Uses XGBoost binary classification (variant vs reference)
- Leave-one-plate-out cross-validation (4 plates per batch)
- Both standard GFP and GFP-intensity-corrected classifications
- Achieves AUROC ~0.976 for both modes

### Data Structure

- **Reference wells**: A19, I21 (CCM2 reference)
- **Variant wells**: M19, C23 (CCM2_Val53Ile)
- **Batches**: 18 and 19 (4 plates each = 8 plates total)
- **Classifiers**: 32 total (16 per batch, 4 well pairs × 4 plates)

## Usage

### Classification

```bash
# Run XGBoost classification
pixi run classify
```

Output files in `data/interim/classification_results/`:
- `classifier_info.csv` - Per-classifier AUROC scores
- `feature_importance.csv` - Feature weights
- `predictions.parquet` - Cell-level predictions

Summary in `data/processed/ccm2_val53ile_phenotype_score.csv`.

### Data Extraction (one-time setup)

If using data from VarChAMP repo:

```bash
# Extract CCM2 profiles from VarChAMP
pixi run extract-data

# Extract cell locations for visualization
pixi run extract-locations
```

### Manuscript Figures

Generate publication-ready figures for the CCM2_Val53Ile analysis:

```bash
# Feature importance barplot (top 10 XGBoost features) -> FigA_feat_importance.png
pixi run feature-importance

# Cell crop figures
pixi run cell-crops          # Reproduce mode: 5 cells/allele, exact notebook cells
pixi run cell-crops-feature  # Feature-based: 50 cells/allele, percentile selection
pixi run cell-crops-all      # Generate both figures
```

**Cell Crop Figure Modes:**

| Mode | Output | Description |
|------|--------|-------------|
| `reproduce` | `FigB_reproduce_cell_crops.png` | 5 cells × 2 alleles, exact cells from reference notebook |
| `feature` | `FigB_feat_sel_cell_crops.png` | 50 cells × 2 alleles, feature-based selection |

Feature-based selection uses `Cytoplasm_Intensity_MinIntensityEdge_GFP`:
- Lower-mean allele: cells from 10-50% percentile
- Higher-mean allele: cells from 50-90% percentile

All cell crop figures include 20 µm scale bars (Phenix 20x, 2×2 binning = 0.598 µm/pixel).

Output directory: `data/processed/manuscript_figures/`

### Standalone Image Download

To download images without VarChAMP repo access:

```bash
# Download all CCM2 well images (default: sites 1-5)
pixi run download-images

# Download specific wells/sites
pixi run download-images -- --wells A19,M19 --sites 1,2,3
```

## Directory Structure

```
.
├── data/
│   ├── raw/
│   │   └── single_cell_profiles/    # Extracted cell profiles (parquet)
│   ├── interim/
│   │   ├── classification_results/  # Classifier outputs
│   │   └── cell_locations/          # Cell coordinates for cropping
│   ├── processed/
│   │   └── manuscript_figures/      # Publication-ready figures
│   └── raw/images/                  # Downloaded TIFF images (optional)
├── scripts/
│   ├── 01_xgb_classification.py     # Main classification script
│   ├── 02_feature_importance_figure.py  # Feature importance barplot
│   ├── 03_cell_crop_figures.py      # Cell crop visualization figures
│   ├── extract_ccm2_profiles.py     # Profile extraction from VarChAMP
│   ├── extract_cell_locations.py    # Location extraction for visualization
│   └── download_cpg_images.py       # CPG image downloader
└── pixi.toml                        # Project configuration
```

## Reproduction Workflow

### From Scratch (requires VarChAMP access)

```bash
pixi run extract-data
pixi run extract-locations
pixi run classify
pixi run feature-importance
pixi run cell-crops-all
```

### Using Provided Data

```bash
pixi run classify
pixi run feature-importance
pixi run cell-crops-all
```

## Results Summary

| Metric | Standard GFP | GFP-Adjusted |
|--------|-------------|--------------|
| AUROC | 0.9760 ± 0.0359 | 0.9801 ± 0.0310 |
| Batch 18 | 0.9710 ± 0.0472 | 0.9762 ± 0.0414 |
| Batch 19 | 0.9810 ± 0.0195 | 0.9840 ± 0.0156 |
| Classifiers | 32 | 32 |

## Dependencies

Managed via pixi (see `pixi.toml`):
- Python 3.10+
- XGBoost 2.0+
- Polars, Pandas, NumPy
- scikit-learn
- tifffile
