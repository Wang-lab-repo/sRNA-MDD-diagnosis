# Plasma Small RNA Classifier for Major Depressive Disorder

[![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](LICENSE)

**Development and Prospective Validation of a Plasma Small RNA Classifier for Major Depressive Disorder**



---

## Overview

This repository contains the complete analysis pipeline for a multicenter diagnostic study (N = 1,124) that developed and prospectively validated a plasma small RNA (sRNA) classifier for distinguishing major depressive disorder (MDD) from non-MDD participants across 5 clinical centers in China. The locked 21-sRNA Random Forest classifier was evaluated in 2 retrospective external and 2 prospectively recruited validation cohorts without refitting.


---

## Repository Structure

```
.
├── README.md
├── LICENSE
├── requirements.txt              # Python dependencies
├── Environment.md                # Full software environment & external tools
├── scripts/
│   ├── 01_preprocessing/         # Small RNA sequencing preprocessing (QC, trimming, annotation)
│   │   ├── sRNA_QC_Pipeline.md
│   │   ├── miRNA_sRNA_Pipeline.md
│   │   ├── SPORTS_sRNA_Pipeline.md
│   │   ├── piRNA_sRNA_Pipeline.md
│   │   ├── filter_by_quality.py
│   │   ├── mirna_merge.py
│   │   └── sports_merge.py
│   ├── 02_data_integration/      # CPM normalization, cohort splitting
│   │   └── integrate_batches.py
│   ├── 03_differential_expression/  # Differential expression analysis
│   │   ├── DE_code.md
│   │   ├── median_mirna.py
│   │   ├── median_pirna.py
│   │   ├── median_rsrna.py
│   │   └── median_tsrna.py
│   ├── 04_model_training/        # Random Forest classifier training & feature selection
│   │   └── train_classifier.py
│   ├── 05_evaluation/            # External validation with Platt scaling & bootstrap CI
│   │   └── evaluate.py
│   ├── 06_figures/               # ROC, DCA, and calibration plots
│   │   └── plot_figures.py
│   └── 07_pathway_analysis/      # Target prediction, GO/KEGG enrichment
│       ├── enrichgo.py
│       ├── intersect_targets.py
│       ├── parse_rnahybrid.py
│       ├── run_miranda_parser.py
│       ├── GO_Enrichment_Bar_Plot.r
│       ├── Functional_Enrichment_Analysis_Workflow.md
│       └── sRNA_Target_Prediction_Pipeline.md
└── References
```

---

## Pipeline Summary

### 1. Small RNA Preprocessing (`scripts/01_preprocessing/`)
- **Quality control**: FastQC v0.12.1 + Trim Galore v0.6.7 (Q >= 25, length >= 18 nt)
- **Advanced filtering**: Base-quality filtering via `filter_by_quality.py`
- **miRNA annotation**: miRDeep2 v2.0.1.3 against miRBase (hg19)
- **tsRNA/rsRNA annotation**: SPORTS1.1 against GtRNAdb, rRNAdb, piRBase
- **piRNA quantification**: Bowtie v1.0.0 + RSEM against piRBase
- After QC: 3,009 features retained (454 miRNAs, 2,170 rsRNAs, 385 tsRNAs)

### 2. Data Integration (`scripts/02_data_integration/`)
- Merge miRNA, rsRNA, tsRNA expression across two sequencing batches
- CPM normalization
- Feature filtering: retain features with >=30% non-zero expression in the derivation cohort
- Split into three analytic datasets: derivation, external validation, prospective validation

### 3. Differential Expression (`scripts/03_differential_expression/`)
- Median-based fold-change calculation
- Two-sided Mann-Whitney U test with Benjamini-Hochberg FDR correction
- Volcano plot visualization

### 4. Classifier Training (`scripts/04_model_training/`)
- **Feature selection** (performed exclusively within the derivation cohort):
  1. Random Forest Gini importance (top 15)
  2. Mean absolute SHAP importance, all derivation-cohort participants (top 15)
  3. SHAP importance restricted to MDD vs. other psychiatric disease comparisons (top 15)
  - Final panel: union of top-15 features from all three rankings -> 21 sRNA features
- **Model**: Random Forest with class-balanced subsampling
- **Hyperparameter tuning**: RandomizedSearchCV with 5-fold stratified CV (200 iterations, AUC scoring)
- All steps locked before any validation data were analyzed

### 5. External Validation (`scripts/05_evaluation/`)
- Platt scaling (logistic calibration) fitted on derivation-cohort probabilities
- Classification threshold (Youden-optimal) locked from derivation cohort
- Per-center evaluation with bootstrap 95% confidence intervals (1,000 iterations)
- Metrics: AUC, AUPRC, Accuracy, Sensitivity, Specificity, PPV, NPV, F1

### 6. Figures (`scripts/06_figures/`)
- ROC curves with AUC
- Decision Curve Analysis (DCA)
- Calibration plots (Brier score, calibration slope, calibration-in-the-large)

### 7. Pathway Analysis (`scripts/07_pathway_analysis/`)
- Target gene prediction (miRanda v3.3a, RNAhybrid v2.1.2)
- GO enrichment analysis (Fisher's exact test)
- KEGG pathway enrichment
- Network visualization (Cytoscape v3.10.3)

---

## Environment Setup

```bash
# Python environment
conda create -n mdd-srna python=3.11
conda activate mdd-srna
pip install -r requirements.txt
```

For R-based analyses and external bioinformatics tools (miRDeep2, SPORTS1.1, Bowtie, miRanda, RNAhybrid, Cytoscape), see [Environment.md](Environment.md).

---

## Usage

Scripts are designed to run sequentially. The training and evaluation scripts expect data prepared by the integration step:

```bash
# Step 1: Data integration (requires batch1/ and batch2/ input directories)
python scripts/02_data_integration/integrate_batches.py

# Step 2: Train classifier
python scripts/04_model_training/train_classifier.py

# Step 3: Evaluate on external validation cohorts
python scripts/05_evaluation/evaluate.py

# Step 4: Generate figures
python scripts/06_figures/plot_figures.py
```

---

## Data Availability

The individual-level sequencing and clinical data that support the findings of this study are available from the corresponding author upon reasonable request, subject to institutional ethics committee approval and data-sharing agreements. The derived feature matrices and model are not publicly deposited due to patient privacy constraints.


---

## License

This project is licensed under the Apache License 2.0. See [LICENSE](LICENSE) for details.
