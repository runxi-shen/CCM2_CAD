# Cell Painting Analysis
Cell painting was used to understand the imapct of CCM2 V53I variant. 

In U2OS cells, CCM2 WT/V53I-GFP were expressed. CCM2 localization features and U2OS morphology features were extracted by CellProfiler for analysis.
Similarly, WT/V53I CCM2-V5 was expressed in teloHAEC cells, we trained a new CellProfiler pipeline to extract the teloHAEC morphology and CCM2 localization features. 

Measurements for each feature under each cell compartments (e.g., V53IEC_Cells.csv, V53IEC_Cytoplasm.csv, V53IEC_Nuclei.csv) are stored under DATA_DIR. Each row represents a single segmented object (e.g., a cell) and columns are morphological features extracted by CellProfiler. Metadata columns indicate experimental group/genotype is encoded in Metadata_Genotype where W prefixes denote WT and V prefixes denote CCM2-V53I.
The full data tables are not included in this repository due to file size limitations. Contact the author for processed data files.

