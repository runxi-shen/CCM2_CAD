# CCM2 Heterozygous Knockout — Analysis Repository

This repository contains analysis code accompanying the manuscript on CCM2 heterozygous knockout and atherosclerosis development, integrating single-cell RNA sequencing and cell painting assays.

## Repository Structure

```
├── scRNAseq/                      ← scRNA-seq analysis (in vivo, mice aorta)
│   ├── scripts/                   ← Preprocessing (CellBender, Scanpy, CellRank)
│   ├── analysis/                  ← Clustering, annotation, EC subclustering (R)
│   ├── data/                      ← Gene sets & data documentation
│   ├── results/                   ← Output tables
│   └── figures/
│
├── cell_painting/                 ← Cell painting analysis (in vitro)
│   ├── teloHAEC/                  ← teloHAEC analysis (this repo)
│   │   ├── pipeline/              ← CellProfiler pipeline files (.cppipe)
│   │   ├── analysis/              ← Data analysis notebooks (Python)
│   │   ├── results/
│   │   └── figures/
│   └── README.md                  ← Overview + link to U2OS analysis
```

## scRNA-seq

Single-cell RNA sequencing of aortic tissue from CCM2 heterozygous knockout and wild-type mice at Week 0 and Week 12. See [`scRNAseq/README.md`](scRNAseq/README.md) for details.

Raw data: [GEO GSE315884](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE315884)

## Cell Painting

Morphological profiling using the Cell Painting assay in two in vitro systems. The teloHAEC analysis is provided in this repository; the U2OS analysis was performed by collaborators and is available separately. See [`cell_painting/README.md`](cell_painting/README.md) for details.

## Contact

Shi Fang — [GitHub](https://github.com/HouseFang-Broad)
