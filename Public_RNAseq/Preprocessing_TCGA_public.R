################################################################################
# TCGA LUAD Data Preparation Pipeline for VAE
# 
# Complete workflow:
#   1. Download & load TCGA-LUAD data
#   2. Extract primary tumors
#   3. Map Ensembl IDs to gene symbols
#   4. Remove duplicate symbols (keep highest MAD)
#   5. Select top 5000 variable genes
#   6. Create condition matrix (smoking status + age group)
#   7. Build VAE input matrices
#   8. Export to CSV
################################################################################

library(TCGAbiolinks)
library(SummarizedExperiment)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(edgeR)
library(matrixStats)
library(tidyverse)

setwd("Documents/digital_patients/")

# ==============================================================================
# 1. Download TCGA-LUAD Data
# ==============================================================================

query_luad <- GDCquery(
  project = "TCGA-LUAD",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

GDCdownload(query_luad, method = "api", files.per.chunk = 10)
tcga_data <- GDCprepare(query_luad)
tcga_data %>% saveRDS("rds/TCGA_LUAD_data.rds")

# ==============================================================================
# 2. Extract Data & Filter to Primary Tumors
# ==============================================================================

tcga_data <- readRDS("rds/TCGA_LUAD_data.rds")

tcga_counts <- assay(tcga_data, "unstranded")
tcga_tpm <- assay(tcga_data, "tpm_unstrand")
tcga_metadata <- colData(tcga_data) %>% as.data.frame()

# Keep only primary tumor samples
tcga_metadata <- tcga_metadata %>%
  filter(sample_type == "Primary Tumor")

primary_barcodes <- rownames(tcga_metadata)
tcga_counts <- tcga_counts[, primary_barcodes]
tcga_tpm <- tcga_tpm[, primary_barcodes]

# ==============================================================================
# 3. Map Ensembl IDs to Gene Symbols
# ==============================================================================

# Remove version numbers from Ensembl IDs
gene_ids_clean <- sub("\\..*", "", rownames(tcga_tpm))

# Map to gene symbols
gene_symbols <- mapIds(org.Hs.eg.db,
                       keys = gene_ids_clean,
                       column = "SYMBOL",
                       keytype = "ENSEMBL",
                       multiVals = "first")

# Keep only genes with symbol mapping
has_symbol <- !is.na(gene_symbols)
tcga_tpm_symbol <- tcga_tpm[has_symbol, ]
tcga_counts_symbol <- tcga_counts[has_symbol, ]
rownames(tcga_tpm_symbol) <- gene_symbols[has_symbol]
rownames(tcga_counts_symbol) <- gene_symbols[has_symbol]

# ==============================================================================
# 4. Handle Duplicate Gene Symbols (keep highest MAD)
# ==============================================================================

dup_symbols <- rownames(tcga_tpm_symbol)[duplicated(rownames(tcga_tpm_symbol)) | 
                                         duplicated(rownames(tcga_tpm_symbol), fromLast = TRUE)]
unique_dup <- unique(dup_symbols)

keep_idx <- !duplicated(rownames(tcga_tpm_symbol))

for(sym in unique_dup) {
  idx <- which(rownames(tcga_tpm_symbol) == sym)
  if(length(idx) > 1) {
    mads <- rowMads(tcga_tpm_symbol[idx, ])
    keep_idx[idx] <- FALSE
    keep_idx[idx[which.max(mads)]] <- TRUE
  }
}

tcga_tpm_final <- tcga_tpm_symbol[keep_idx, ]
tcga_counts_final <- tcga_counts_symbol[keep_idx, ]

# ==============================================================================
# 5. Select Top 5000 Variable Genes
# ==============================================================================

gene_variances <- rowVars(tcga_tpm_final)
top_n_genes <- 5000
variance_threshold <- sort(gene_variances, decreasing = TRUE)[top_n_genes]
keep_high_var <- gene_variances >= variance_threshold

tcga_tpm_final_hv <- tcga_tpm_final[keep_high_var, ]
tcga_counts_final_hv <- tcga_counts_final[keep_high_var, ]

cat("Primary tumor samples:", nrow(tcga_metadata), "\n")
cat("Genes after symbol mapping:", nrow(tcga_tpm_final), "\n")
cat("Final genes (top 5000):", nrow(tcga_tpm_final_hv), "\n")

# ==============================================================================
# 6. Create Condition Matrix (Smoking Status + Age Group)
# ==============================================================================

# Extract relevant metadata for condition
df <- tcga_metadata %>%
  select(gender, age_at_index, tobacco_smoking_status) %>%
  mutate(
    smoking_status = case_when(
      tobacco_smoking_status == "Lifelong Non-Smoker" ~ "Never",
      grepl("Reformed", tobacco_smoking_status) ~ "Former",
      tobacco_smoking_status == "Current Smoker" ~ "Current",
      TRUE ~ NA_character_
    ),
    age_group = ifelse(age_at_index < 65, "<65", ">=65")
  )

# Filter out samples with missing values
df_filtered <- df %>%
  filter(!is.na(smoking_status) & !is.na(age_group)) %>%
  select(smoking_status, age_group)

cat("Samples after removing missing:", nrow(df_filtered), "\n")

# Create condition matrix with one-hot encoding
meta_cond <- df_filtered %>%
  mutate(
    smoking_status = factor(smoking_status, levels = c("Never", "Former", "Current")),
    age_group = factor(age_group, levels = c("<65", ">=65"))
  )

cond_mat <- model.matrix(~ smoking_status + age_group - 1, data = meta_cond)

# ==============================================================================
# 7. Match Expression & Condition Data, Build VAE Input
# ==============================================================================

sample_order <- rownames(meta_cond)

# Reorder expression data to match condition data
tcga_tpm_final_hv_ordered <- tcga_tpm_final_hv[, sample_order]

# Log transformation & z-score normalization
expr_mat <- log2(tcga_tpm_final_hv_ordered + 1)
expr_mat_scaled <- t(scale(t(expr_mat)))

# Verify alignment
stopifnot(all(colnames(expr_mat_scaled) == rownames(meta_cond)))
stopifnot(nrow(meta_cond) == nrow(cond_mat))

# Prepare VAE input
tcga_vae_input <- list(
  expr = t(expr_mat_scaled),  # samples x genes
  cond = cond_mat,             # samples x condition features
  meta = df_filtered
)

saveRDS(tcga_vae_input, file = "rds/01.tcga_luad_vae_input.rds")

# ==============================================================================
# 8. Export to CSV
# ==============================================================================

expr_df <- as.data.frame(tcga_vae_input$expr)
cond_df <- as.data.frame(tcga_vae_input$cond)

write.csv(expr_df, "rds/tcga_luad_expr.csv", quote = FALSE)
write.csv(cond_df, "rds/tcga_luad_cond.csv", quote = FALSE)

cat("\n✓ VAE input saved to RDS and CSV\n")
