# Data

## Raw Data Availability

Raw sequencing data (FASTQ files) and processed CellRanger outputs are available from GEO:

**GEO Accession:** [GSE315884](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE315884)

## Input Files for R Analysis

The R markdown scripts expect the following preprocessed files in this directory:

### Preprocessed Seurat Objects (from Python pipeline)

These RDS files are converted from h5ad (AnnData) objects after CellBender ambient RNA removal and Scanpy QC:

| Filename | Sample | Genotype | Timepoint |
|----------|--------|----------|-----------|
| `CCM2WTWeek0_preprocessed.rds` | CCM2WTWeek0 | Wild-type | Week 0 (Baseline) |
| `CCM2HetWeek0_preprocessed.rds` | CCM2HetWeek0 | CCM2 Het KO | Week 0 (Baseline) |
| `CCM2WTWeek12_preprocessed.rds` | CCM2WTWeek12 | Wild-type | Week 12 (PCSK9) |
| `CCM2HetWeek12_preprocessed.rds` | CCM2HetWeek12 | CCM2 Het KO | Week 12 (PCSK9) |

### CellRanger H5 Files (available from GEO)

Raw CellRanger filtered feature-barcode matrices:

- `CCM2WTWeek0_filtered_feature_bc_matrix.h5`
- `CCM2HetWeek0_filtered_feature_bc_matrix.h5`
- `CCM2WTWeek12_filtered_feature_bc_matrix.h5`
- `CCM2HetWeek12_filtered_feature_bc_matrix.h5`

### Custom Gene Sets

- `Custom_gene_sets_Mm.csv` - Custom gene sets for module score analysis (included in this repository)

## Preprocessing Pipeline

The Python preprocessing pipeline (in `scripts/`) performs:

1. **CellBender** - Ambient RNA removal
2. **Scanpy QC** - Including:
   - Scrublet doublet detection
   - Cell filtering: nUMI > 500, nGene > 250, log10GenesPerUMI > 0.8, mitoRatio < 0.2
3. **Conversion** - h5ad to RDS for Seurat analysis

## Experimental Design

- **Species:** Mus musculus (mouse)
- **Tissue:** Aorta
- **Model:** CCM2 heterozygous knockout (Het) vs Wild-type (WT)
- **Timepoints:** 
  - Week 0: Baseline
  - Week 12: PCSK9-induced atherosclerosis
- **Sequencing:** 10x Genomics 3' scRNA-seq
