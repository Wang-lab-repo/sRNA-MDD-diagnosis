# Environment.md
**Project**: *sRNA-MDD-Diagnosis*  
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

| Package        | Description                                     |
| -------------- | ----------------------------------------------- |
| `numpy`        | Numerical computing                             |
| `pandas`       | Data wrangling and I/O                          |
| `scikit-learn` | Machine learning algorithms and evaluation      |
| `xgboost`      | Gradient boosting classifier                    |
| `lightgbm`     | Gradient boosting framework                     |
| `catboost`     | Gradient boosting with categorical features     |
| `shap`         | SHAP values for model explainability            |
| `matplotlib`   | Plotting and visualization                      |
| `seaborn`      | Statistical data visualization                  |
| `pymoo`        | Multi-objective optimization toolkit (NSGA-III) |
| `joblib`       | Model serialization and parallel processing     |
| `scipy`        | Scientific computing and statistical functions  |
| `statsmodels`  | Statistical modeling and inference              |

---

## 2. Jupyter Notebook Environment

Jupyter notebooks are run within the configured Python environment.

```bash
jupyter lab
```

Ensure the following package is installed: `notebook` or `jupyterlab`.

---

## 3. R Environment

- **Version**: R 4.4.1 (recommended)

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

## 4. Small RNA Annotation & Preprocessing Tools

### miRDeep2
- **Version**: v2.0.1.3
- **Purpose**: Identification of miRNAs (precursor and mature)
- **Platform**: Linux command-line (Perl required)

### Bowtie
- **Version**: 1.0.0
- **Purpose**: Sequence alignment for piRNA (tsRNA, rsRNA)
- **Platform**: Linux command-line, integrated in SPORTS1.1

### SPORTS1.1
- **Purpose**: Annotation of tsRNA and rsRNA
- **Dependencies**: Bowtie, Perl, Python 3

---

## 5. Target Prediction & Enrichment Analysis

### miRanda
- **Version**: v3.3a (based on Python 3.9)
- **Purpose**: miRNA target prediction

### RNAhybrid
- **Version**: v2.1.2 (based on Python 3.9)
- **Purpose**: Free energy-based target binding prediction

---

## 6. Network Visualization

### Cytoscape
- **Version**: v3.10.3
- **Purpose**: Visualization of gene-sRNA-target-GO networks
- **Platform**: Desktop GUI for Windows/macOS/Linux
