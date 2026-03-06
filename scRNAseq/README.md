# CCM2hetKOmice_aorta_scRNAseq
Aortas CCM2het mice at baseline or with atherosclerosis
# Data Directory

## Data Availability

Raw and processed data are available from GEO: [GSE315884](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE315884)

## Samples

| Sample | Genotype | Timepoint | Description |
|--------|----------|-----------|-------------|
| CCM2WTWeek0 | Wild-type | Week 0 | Baseline WT |
| CCM2HetWeek0 | CCM2 Het KO | Week 0 | Baseline Het |
| CCM2WTWeek12 | Wild-type | Week 12 | Atherosclerosis WT |
| CCM2HetWeek12 | CCM2 Het KO | Week 12 | Atherosclerosis Het |

## Preprocessing Pipeline

### Step 1: CellBender (Ambient RNA Removal)
- **Script:** `scripts/01_cellbender.ipynb`
- **Input:** CellRanger `filtered_feature_bc_matrix.h5`
- **Output:** `*_cellbender_output_filtered.h5`

### Step 2: Scanpy QC + Scrublet
- **Script:** `scripts/02_scanpy_qc_scrublet.ipynb`
- **Input:** CellBender output
- **Output:** `*_qc_scrub.h5ad`

### Step 3: Convert to Seurat
- **Script:** `scripts/03_h5ad_to_seurat.R`
- **Input:** QC'd h5ad files
- **Output:** `*_preprocessed.rds`

## Files Required for Analysis

Place these files in the `data/` directory:

### For `01_clustering_annotation.Rmd`
- `CCM2WTWeek0_preprocessed.rds`
- `CCM2HetWeek0_preprocessed.rds`
- `CCM2WTWeek12_preprocessed.rds`
- `CCM2HetWeek12_preprocessed.rds`

### For `02_EC_subclustering.Rmd`
- `Custom_gene_sets_Mm.csv` (custom gene sets for module scoring)

### For `04_cellrank_trajectory.ipynb`
- `EC_annotated.h5ad` (exported from Seurat EC object)

### For `05_cellrank_sankey.Rmd`
- `cellrank_fate_counts.xlsx` (output from CellRank analysis)

## Output Files Generated

### From `01_clustering_annotation.Rmd`
- `results/filterd_all_annotated.qs` - All cells, annotated
- `results/ECs_processed.qs` - EC subset, recentered

### From `02_EC_subclustering.Rmd`
- `results/EC_annotated.qs` - EC with subcluster annotations
- `results/EC_analysis/longnumCells_EC.csv` - Subcluster counts
- `results/EC_analysis/customized_setscore_ECs_*.csv` - Gene set scores

