# ==============================================================================
# 3. Convert h5ad to Seurat RDS
# ==============================================================================
#
# This script converts QC'd h5ad files (from Scanpy/Scrublet) to Seurat RDS format.
#
# Input: *_qc_scrub.h5ad files (from 02_scanpy_qc_scrublet.ipynb)
# Output: *_preprocessed.rds files (for R analysis)
#
# ==============================================================================

# Load libraries
library(anndata)
library(Seurat)
library(SeuratDisk)

# ------------------------------------------------------------------------------
# Set paths - UPDATE THESE FOR YOUR ENVIRONMENT
# ------------------------------------------------------------------------------

# Directory containing h5ad files
input_dir <- "./"

# Directory for output RDS files (can be same as input)
output_dir <- "./"

# ------------------------------------------------------------------------------
# Sample names
# ------------------------------------------------------------------------------

samples <- c(
  "CCM2WTWeek0",
  "CCM2HetWeek0",
  "CCM2WTWeek12",
  "CCM2HetWeek12"
)

# ------------------------------------------------------------------------------
# Convert h5ad to Seurat RDS
# ------------------------------------------------------------------------------

for (sample in samples) {
  
  cat("\n", paste(rep("=", 60), collapse = ""), "\n")
  cat("Processing:", sample, "\n")
  cat(paste(rep("=", 60), collapse = ""), "\n")
  
  # Define file paths
  h5ad_path <- file.path(input_dir, paste0(sample, "_qc_scrub.h5ad"))
  rds_path <- file.path(output_dir, paste0(sample, "_preprocessed.rds"))
  
  # Check if input file exists
  if (!file.exists(h5ad_path)) {
    cat("  ERROR: Input file not found:", h5ad_path, "\n")
    next
  }
  
  # Read h5ad file
  cat("  Reading h5ad file...\n")
  h5ad_obj <- read_h5ad(h5ad_path)
  
  # Create Seurat object
  # Note: h5ad stores data as cells x genes, Seurat expects genes x cells
  cat("  Creating Seurat object...\n")
  seurat_obj <- CreateSeuratObject(
    counts = t(as.matrix(h5ad_obj$X)),
    meta.data = h5ad_obj$obs
  )
  
  # Report cell/gene counts
  cat("  Cells:", ncol(seurat_obj), "\n")
  cat("  Genes:", nrow(seurat_obj), "\n")
  
  # Report metadata columns transferred
  cat("  Metadata columns:", paste(colnames(seurat_obj@meta.data), collapse = ", "), "\n")
  
  # Save as RDS
  cat("  Saving RDS file...\n")
  saveRDS(seurat_obj, file = rds_path)
  
  cat("  Saved to:", rds_path, "\n")
}

cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("Conversion complete!\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

# ------------------------------------------------------------------------------
# Output files ready for R analysis:
# - CCM2WTWeek0_preprocessed.rds
# - CCM2HetWeek0_preprocessed.rds
# - CCM2WTWeek12_preprocessed.rds
# - CCM2HetWeek12_preprocessed.rds
#
# These files are loaded by analysis/01_clustering_annotation.Rmd
# ------------------------------------------------------------------------------
