# [Research Proposal] Development of Lung Cancer-Specific Deconvolution Methodology and TME-Index for Prognostic Stratification

## 1. Background and Rationale

### 1.1 Limitations of Current TNM Staging System

The TNM staging system, the standard diagnostic framework for lung cancer, relies on anatomical features and fails to capture the biological heterogeneity underlying patient-to-patient prognostic variability. Patients within the same stage often exhibit drastically different recurrence patterns, highlighting the need for objective molecular biomarkers that complement morphological staging.

### 1.2 Observed Heterogeneity in Transcriptomic Data Distribution

Unsupervised dimensionality reduction (UMAP) analysis of TCGA-LUAD data reveals that patient samples do not cluster distinctly by clinical stage (Stage 1-4) but rather display broad overlap across stages. This suggests that global gene expression patterns alone cannot fully explain the existing staging system, and that extraction of key variables such as tumor microenvironment (TME) composition is necessary for precise analysis.

### 1.3 Limitations of Existing General-Purpose Deconvolution Tools

Existing general-purpose deconvolution tools such as CIBERSORTx, EPIC, and MuSiC utilize reference profiles derived from blood cells or pan-cancer datasets, which inadequately capture lung cancer-specific stromal and immune suppressive environments. Notably, lung cancer exhibits critical heterogeneity in cancer-associated fibroblast (CAF) subtypes (inflammatory CAF, myofibroblastic CAF) and immunosuppressive myeloid populations that significantly impact prognosis, yet these nuances are poorly resolved by generic references. Furthermore, the diverse algorithmic approaches (SVR, NNLS, weighted least squares) employed by these tools necessitate systematic comparison and validation for optimal application to lung cancer transcriptomic data.

---

## 2. Research Objectives and Methods

### 2.1 Research Objectives

* **Development of lung cancer-specific deconvolution methodology:** Establishment of lung cancer-tailored single-cell references and optimization of existing algorithms, with potential ensemble approaches, through systematic algorithm comparison
* **Derivation of TME-Index:** High-resolution quantification of cancer, immune, and stromal cell proportions using the developed methodology
* **Prognostic biomarker research:** Analysis of discordance between staging and biological state, integrated with survival data to establish TME-Index as a complementary prognostic indicator

### 2.2 Research Methods

#### [Phase 1-1] Construction of Lung Cancer-Specific Single-Cell Reference

Extract core gene signatures for cancer, immune, and stromal cell types directly from **NSCLC-specific single-cell data in TISCH (Tumor Immune Single-cell Hub)**. Unlike generic references, this approach constructs a high-precision signature matrix that reflects lung cancer microenvironment-specific cell types and states.

#### [Phase 1-2] Systematic Algorithm Comparison and Methodology Optimization

Systematically compare the performance of multiple deconvolution algorithms using the constructed lung cancer-specific reference:

- **QP (Quadratic Programming):** Explicitly enforces the physical constraint that cell proportions sum to 100%
- **NNLS (Non-Negative Least Squares):** Constrained least squares approach preventing negative proportions
- **SVR (Support Vector Regression):** Core algorithm of CIBERSORTx, nu-regression based
- **Robust Regression:** Regression methods robust to outliers

**Performance Evaluation Metrics:**
- **Accuracy:** MAE (Mean Absolute Error), RMSE, CCC (Concordance Correlation Coefficient)
- **Biological validity:** Correlation between marker gene expression and estimated cell proportions
- **Reconstruction accuracy:** Fidelity of pseudobulk reconstruction to original expression patterns
- **Prognostic power:** C-index based on TCGA survival data

Comprehensively evaluate accuracy, stability, and computational efficiency of each algorithm, and establish a lung cancer-optimized deconvolution pipeline through ensemble or constraint optimization as needed. The goal is not to develop an entirely novel algorithm, but rather to readjust existing algorithms to lung cancer-specific references and constraints, identifying the optimal combination.

#### [Phase 2] Stage-Specific Distribution Analysis and Outlier Identification

Apply the developed methodology to the TCGA cohort to establish standard ranges of Cancer:Immune:Stroma ratios for each stage. Through statistical analysis, identify 'stage-state discordant' patient subgroups who exhibit low pathological stage but high-risk biological patterns.

#### [Phase 3] Validation of Prognostic Power Using Survival Data

Integrate the derived TME-Index with actual TCGA survival data (Overall Survival). Compare how well the existing staging system versus the TME-Index explains patient prognosis through **Kaplan-Meier analysis and C-index**.

---

## 3. Project Roadmap

| Phase | Main Research Activities | Expected Deliverables |
|-------|-------------------------|----------------------|
| **Phase 1-1** | Construction of TISCH-based lung cancer-specific single-cell signature matrix | Lung cancer-specific reference signatures |
| **Phase 1-2** | Performance comparison of deconvolution algorithms (QP, NNLS, SVR, etc.) and methodology optimization | Lung cancer-optimized deconvolution pipeline |
| **Phase 2** | Stage-specific TME distribution modeling and discordant patient identification | Statistical analysis model of stage-TME discordance |
| **Phase 3** | Validation of prognostic power using survival data | Index validation and comparative data |
| **Phase 4** | Indicator development and prognostic scoring logic establishment | Prognostic scoring prototype |

---

## 4. Research Strengths and Expected Impact

### 4.1 Enhanced Reliability Through Systematic Methodological Comparison

Unlike existing studies that rely on single tools or algorithms, this research directly compares and validates multiple deconvolution algorithms applied to lung cancer-specific references using quantitative metrics (MAE, CCC, marker correlation, C-index). This establishes a lung cancer-optimized pipeline while ensuring objectivity and reproducibility of results. Rather than developing an entirely new algorithm, this study aims to readjust existing methodologies to lung cancer-specific conditions and identify optimal combinations.

### 4.2 High Analytical Consistency and Feasibility

This research ensures analytical consistency by utilizing lung cancer-specific single-cell data rather than generic references. By employing **core cell proportions** as intuitive variables instead of complex black-box modeling, results are readily interpretable and the approach offers high feasibility for immediate integration into existing transcriptomic analysis environments.

### 4.3 Value as Clinical Decision Support Indicator

By providing data-driven evidence of the gap between morphological staging and actual biological prognosis, this research establishes groundwork for future clinical application as a complementary prognostic indicator. This particularly contributes to developing personalized management strategies for identifying high-risk recurrence groups among early-stage lung cancer patients.

---

## Addendum: Anticipated Questions and Answers

### Q1. Why develop a new methodology instead of using existing tools like CIBERSORTx or EPIC?

**A:** Existing tools use generic references (blood cells or pan-cancer) that fail to precisely capture lung cancer-specific microenvironments. In particular, they struggle to distinguish critical lung cancer features such as CAF subtypes (inflammatory/myofibroblastic) or heterogeneity in immunosuppressive myeloid populations. This research constructs references from TISCH lung cancer-specific single-cell data and quantitatively compares multiple algorithms to establish a pipeline optimized for lung cancer data.

### Q2. What defines an "optimal methodology"?

**A:** We comprehensively evaluate using four metrics:
- **Accuracy:** Prediction error measured by MAE, RMSE, CCC
- **Biological validity:** Correlation between marker gene expression and estimated cell proportions
- **Reconstruction accuracy:** Fidelity of pseudobulk reconstruction to original transcriptome
- **Prognostic power:** C-index based on TCGA survival data

These metrics identify algorithmic strengths and weaknesses, enabling selection of the most stable and accurate methodology for lung cancer data.

### Q3. Are you developing a completely new algorithm?

**A:** No. This research aims to readjust existing algorithms (QP, NNLS, SVR, Robust Regression) to lung cancer-specific references and constraints, potentially deriving optimal combinations through ensemble approaches. The focus is on adapting and optimizing validated methodologies for lung cancer data rather than developing entirely new mathematical frameworks.

### Q4. Is validation with TCGA data alone sufficient?

**A:** TCGA-LUAD represents a large-scale cohort suitable for methodology development and primary validation. However, this research scope focuses on methodology establishment and proof-of-concept, recognizing that additional validation in independent cohorts (GEO, ICGC) will be necessary. At the current stage, we prevent overfitting by splitting TCGA into training/test sets internally.

### Q5. Can TME-Index replace the existing staging system?

**A:** The goal is **complementation**, not replacement. The TNM staging system remains the standard providing anatomical information, while TME-Index adds molecular biological information to explain prognostic heterogeneity within stages. In clinical practice, combining both information sources is expected to enable more precise risk stratification.

### Q6. What are the computational costs and feasibility?

**A:** Deconvolution employs lightweight linear algebra-based algorithms executable within minutes on standard transcriptomic analysis pipelines. Compared to deep learning models, hardware requirements are minimal, and results are intuitively interpretable (three cell proportion variables), presenting low barriers to clinical adoption. The approach is immediately integrable into existing RNA-seq analysis environments via R/Python scripts.

### Q7. Does this include immunotherapy response prediction?

**A:** The primary goal is TME-Index establishment and prognostic power validation. Immunotherapy response prediction is set as an extension in Phase 4, requiring additional analysis linking TME-Index (particularly Immune proportions) with treatment response data. The current stage focuses on identifying stage-TME discordant patients and survival stratification.