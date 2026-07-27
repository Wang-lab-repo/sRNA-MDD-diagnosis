# Environment.md
**Project**: *Development and Prospective Validation of a Plasma Small RNA Classifier for Major Depressive Disorder*
**Purpose**: Document software tools, environments, and dependencies used in the analysis pipeline for reproducibility.

---

## 1. Python Environment

- **Version**: Python 3.11
- **Recommended tool**: conda or venv

```bash
conda create -n mdd-srna python=3.11
conda activate mdd-srna
pip install -r requirements.txt
```

### Key Python Packages

| Package          | Description                                          |
| ---------------- | ---------------------------------------------------- |
| `numpy`          | Numerical computing                                  |
| `pandas`         | Data wrangling and I/O                               |
| `scikit-learn`   | Random Forest classifier, Platt scaling, metrics     |
| `shap`           | SHAP values for model explainability                 |
| `matplotlib`     | Plotting and visualization                           |
| `seaborn`        | Statistical data visualization                       |
| `scipy`          | Scientific computing and statistical functions       |
| `statsmodels`    | Statistical modeling and inference                   |
| `joblib`         | Model serialization and parallel processing          |

---

## 2. R Environment

- **Version**: R >= 4.4.1 (recommended)

### Key R Packages

| Package                             | Description                             |
| ----------------------------------- | --------------------------------------- |
| `TxDb.Hsapiens.UCSC.hg38.knownGene` | Gene annotation database                |
| `org.Hs.eg.db`                      | Gene ID mapping                         |
| `clusterProfiler`                   | Gene enrichment analysis                |
| `enrichplot`                        | Enrichment plot visualization           |
| `ggplot2`                           | Data visualization                      |
| `pheatmap`                          | Heatmap visualization                   |

---

## 3. Small RNA Annotation & Preprocessing Tools

### miRDeep2
- **Version**: v2.0.1.3
- **Purpose**: Identification and quantification of known miRNAs
- **Reference**: Friedlander MR, et al. miRDeep2 accurately identifies known and hundreds of novel microRNA genes in seven animal clades. *Nucleic Acids Res*. 2012;40(1):37-52.
- **Platform**: Linux command-line (Perl required)

### Bowtie
- **Version**: 1.0.0
- **Purpose**: Sequence alignment for piRNA quantification
- **Platform**: Linux command-line, integrated in SPORTS1.1

### SPORTS1.1
- **Purpose**: Annotation and quantification of tsRNA and rsRNA
- **Dependencies**: Bowtie, Perl, Python 3

---

## 4. Target Prediction & Enrichment Analysis

### miRanda
- **Version**: v3.3a
- **Purpose**: miRNA target gene prediction

### RNAhybrid
- **Version**: v2.1.2
- **Purpose**: Free energy-based sRNA-target binding prediction

