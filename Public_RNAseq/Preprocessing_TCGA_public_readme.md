# TCGA LUAD Data Preparation for VAE

## Overview

This script provides a complete data preparation pipeline for training a Conditional Variational Autoencoder (CVAE) using TCGA-LUAD (The Cancer Genome Atlas - Lung Adenocarcinoma) gene expression data.

The pipeline downloads raw RNA-seq counts, processes them into expression matrices, and creates condition variables for VAE training.

## Pipeline Steps

### 1. Download TCGA-LUAD Data

Downloads raw RNA-seq data from GDC (Genomic Data Commons) using TCGAbiolinks package.

Input:
- TCGA project: TCGA-LUAD
- Data category: Transcriptome Profiling
- Data type: Gene Expression Quantification (STAR - Counts)

Output:
- rds/TCGA_LUAD_data.rds (raw SummarizedExperiment object)

### 2. Extract Data & Filter to Primary Tumors

Extracts count and TPM (Transcripts Per Million) matrices from the SummarizedExperiment object.

Filters samples to include only primary tumor samples (sample_type == "Primary Tumor").

Processing:
- Extract unstranded counts
- Extract TPM-normalized values
- Filter metadata to primary tumors only
- Subset both expression matrices by primary tumor barcodes

Output:
- tcga_counts: count matrix (genes x samples)
- tcga_tpm: TPM matrix (genes x samples)
- tcga_metadata: filtered metadata

### 3. Map Ensembl IDs to Gene Symbols

Converts Ensembl gene identifiers to human-readable gene symbols.

Steps:
- Remove version numbers from Ensembl IDs (e.g., ENSG00000000003.15 -> ENSG00000000003)
- Use org.Hs.eg.db package to map Ensembl IDs to official gene symbols
- Keep only genes with successful symbol mapping
- Rename row names to gene symbols

Output:
- tcga_tpm_symbol: TPM matrix with gene symbols as row names
- tcga_counts_symbol: counts matrix with gene symbols as row names

### 4. Handle Duplicate Gene Symbols

Multiple Ensembl IDs may map to the same gene symbol. This step resolves duplicates by retaining the gene with the highest Median Absolute Deviation (MAD).

Steps:
- Identify all gene symbols that appear more than once
- For each duplicate symbol, calculate MAD across samples
- Keep the gene with highest MAD (most variable)
- Remove all other copies of the duplicate symbol

Output:
- tcga_tpm_final: deduplicated TPM matrix (genes x samples)
- tcga_counts_final: deduplicated counts matrix (genes x samples)

### 5. Select Top 5000 Variable Genes

Reduces dimensionality by selecting genes with highest variance across all samples.

Steps:
- Calculate variance for each gene across all samples
- Sort genes by variance (descending)
- Select top 5000 genes
- Subset both TPM and counts matrices

Rationale:
- Highly variable genes capture biological signal
- Reduces noise and computational burden
- Standard practice in single-cell and bulk RNA-seq analysis
- 5000 genes is typical for deep learning applications

Output:
- tcga_tpm_final_hv: top 5000 variable genes (5000 x samples)
- tcga_counts_final_hv: top 5000 variable genes (5000 x samples)

### 6. Create Condition Matrix

Builds one-hot encoded condition variables for CVAE training.

Condition variables:
1. Smoking status (3 categories)
   - Never: Lifelong Non-Smoker
   - Former: Reformed smokers
   - Current: Current Smoker

2. Age group (2 categories)
   - <65: age_at_index < 65
   - >=65: age_at_index >= 65

Processing:
- Extract gender, age_at_index, and tobacco_smoking_status from metadata
- Simplify smoking status categories
- Create age groups
- Remove samples with missing values in either condition variable
- One-hot encode both variables using model.matrix()

Output:
- cond_mat: condition matrix (samples x condition_features)
  - Contains one-hot encoded smoking status (3 columns)
  - Contains one-hot encoded age group (2 columns)
  - Total: 5 columns (one column is redundant but retained for clarity)

### 7. Match Expression & Condition Data, Build VAE Input

Aligns expression data with condition matrix and prepares final VAE input.

Steps:
- Extract sample order from condition matrix row names
- Reorder TPM matrix columns to match condition matrix
- Log transform TPM values: log2(TPM + 1)
- Apply gene-wise z-score normalization (zero mean, unit variance)
- Verify sample alignment between expression and condition matrices
- Transpose expression matrix to samples x genes format

Quality checks:
- Verify all column names in expression match row names in condition matrix
- Verify number of samples in condition matrix equals condition feature matrix

Output:
- tcga_vae_input: list containing
  - expr: expression matrix (samples x genes)
  - cond: condition matrix (samples x condition_features)
  - meta: filtered metadata

### 8. Export to CSV

Saves processed data as CSV files for downstream analysis.

Output files:
- rds/tcga_luad_expr.csv: expression matrix (samples x genes)
- rds/tcga_luad_cond.csv: condition matrix (samples x condition_features)
- rds/01.tcga_luad_vae_input.rds: complete RDS object containing all components

## Input Data

The pipeline requires internet access to download data from GDC. No pre-downloaded files are needed.

Working directory: Documents/digital_patients/

## Output Data Specifications

Final VAE input dimensions:
- Expression matrix: (num_samples x 5000)
  - Samples: Primary tumor samples from TCGA-LUAD with complete metadata
  - Genes: Top 5000 genes selected by variance
  - Values: log2-transformed, gene-wise z-score normalized

- Condition matrix: (num_samples x 5)
  - Smoking_statusCurrent: binary (1 if current smoker, 0 otherwise)
  - Smoking_statusFormer: binary (1 if former smoker, 0 otherwise)
  - Smoking_statusNever: binary (1 if never smoker, 0 otherwise)
  - Age_group>=65: binary (1 if age >= 65, 0 otherwise)
  - Age_group<65: binary (1 if age < 65, 0 otherwise)

## Dependencies

Required R packages:
- TCGAbiolinks: Download TCGA data from GDC
- SummarizedExperiment: Handle expression data structures
- org.Hs.eg.db: Gene annotation database
- AnnotationDbi: Gene ID mapping
- edgeR: Normalization (TMM)
- matrixStats: Fast matrix operations
- tidyverse: Data manipulation

Installation:
```r
# Bioconductor packages
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("TCGAbiolinks", "SummarizedExperiment", 
                       "org.Hs.eg.db", "AnnotationDbi", "edgeR", "matrixStats"))

# CRAN packages
install.packages("tidyverse")
```

## Runtime Notes

- First run will download approximately 10-20 GB of data from GDC (may take 1-2 hours)
- Subsequent runs can skip the download step and load from saved RDS files
- Gene symbol mapping may result in loss of 1-5% of genes without valid symbols
- Duplicate symbol resolution typically affects 2-5% of genes
- Final sample count after removing missing values: typically 500-540 samples

## Usage

Run the complete pipeline:
```r
source("01.Data_preparation_LUAD.v1.R")
```

Or load pre-prepared data:
```r
tcga_vae_input <- readRDS("rds/01.tcga_luad_vae_input.rds")
expr <- tcga_vae_input$expr
cond <- tcga_vae_input$cond
meta <- tcga_vae_input$meta
```

## Quality Control

The script includes automatic quality checks:
- Verifies all samples in expression matrix have corresponding condition values
- Confirms sample order consistency across all matrices
- Logs sample counts at each filtering step
- Uses stopifnot() to halt execution if alignment checks fail

## Normalization Details

Gene-wise z-score normalization is applied to prepare data for VAE:

For each gene:
- Calculate mean expression across all samples
- Calculate standard deviation across all samples
- Transform each sample value: z = (value - mean) / std

This normalization:
- Centers each gene at zero
- Scales each gene to unit variance
- Makes genes comparable on the same scale
- Improves VAE training stability

## Data Integrity

Primary data sources:
- TCGA-LUAD project: 540 primary tumor samples
- Gene annotations: NCBI Ensembl (org.Hs.eg.db)
- Expression values: STAR-aligned RNA-seq counts

No manual curation or filtering by external criteria. All filtering based on:
- Sample type (primary tumors only)
- Gene symbol availability (mapped IDs only)
- Gene variability (top 5000 by variance)
- Metadata completeness (no missing condition values)

## References

TCGA-LUAD Project:
https://portal.gdc.cancer.gov/projects/TCGA-LUAD

TCGAbiolinks Package:
https://bioconductor.org/packages/TCGAbiolinks

Gene Annotation:
https://bioconductor.org/packages/org.Hs.eg.db

## Author Notes

This pipeline was designed to prepare TCGA-LUAD data for training conditional VAE models. The choice of 5000 genes and specific condition variables (smoking status, age) can be modified by changing the respective parameters in the script.

For different TCGA projects, modify the project name in Step 1 (e.g., "TCGA-LUSC" for lung squamous cell carcinoma).
