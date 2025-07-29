# 🧪 Environment.md
**Project**: *sRNA-MDD-Diagnosis*  
**Purpose**: Document software tools, environments, and dependencies used in the analysis pipeline for reproducibility.

---

## 1. Python Environment

- **Version**: Python 3.11.5
- **Recommended tool**: conda or venv

```bash
conda create -n mdd-srna python=3.11.5
conda activate mdd-srna
pip install -r requirements.txt
```

### 📦 Key Python Packages

| Package         | Recommended Version | Description                                      |
|----------------|---------------------|--------------------------------------------------|
| `pandas`       | ≥1.5.3              | Data wrangling and I/O                           |
| `numpy`        | ≥1.23               | Numerical operations                             |
| `scikit-learn` | ≥1.0.2              | Machine learning modeling and validation         |
| `xgboost`      | ≥1.7                | XGBoost classifier                               |
| `lightgbm`     | ≥3.3                | LightGBM classifier                              |
| `catboost`     | ≥1.2                | CatBoost classifier                              |
| `shap`         | ≥0.42               | SHAP value-based feature interpretation          |
| `matplotlib`   | ≥3.6                | Plotting and visualization                       |
| `seaborn`      | ≥0.12               | Statistical plots and heatmaps                   |
| `pymoo`        | ≥0.6.1.3            | NSGA-III multi-objective optimization            |
| `scipy`        | ≥1.10               | Statistical testing                              |
| `statsmodels`  | ≥0.13               | Statistical modeling and inference               |

---

## 📒 2. Jupyter Notebook Environment

Jupyter notebooks are run within the configured Python environment.

```bash
jupyter lab
```

Ensure the following package is installed: `notebook` or `jupyterlab`.

---

## 🧬 3. R Environment

- **Version**: R 4.4.1 (recommended)

### 📦 Key R Packages

| Package              | Description                                    |
|----------------------|------------------------------------------------|
| `sva`                | Batch effect correction                        |
| `ggplot2`            | Data visualization                             |
| `clusterProfiler`    | GO/KEGG enrichment analysis                    |
| `org.Hs.eg.db`       | Gene annotation database                       |

---

## 🧬 4. Small RNA Annotation & Preprocessing Tools

### miRDeep2
- **Version**: v2.0.1.2
- **Purpose**: Identification of miRNAs (precursor and mature)
- **Platform**: Linux command-line (Perl required)

### Bowtie
- **Version**: v1.3
- **Purpose**: Sequence alignment for piRNA, tsRNA, rsRNA
- **Platform**: CLI, integrated in SPORTS1.1

### SPORTS1.1
- **Purpose**: Annotation of tsRNA and rsRNA
- **Dependencies**: Bowtie, Perl, Python 2/3

---

## 🔬 5. Target Prediction & Enrichment Analysis

### miRanda
- **Version**: v3.3a
- **Purpose**: miRNA target prediction
- **Platform**: Linux command-line

### RNAhybrid
- **Version**: v2.1.2
- **Purpose**: Free energy-based target binding prediction
- **Platform**: CLI

---

## 🌐 6. Network Visualization

### Cytoscape
- **Version**: v3.10.3
- **Purpose**: Visualization of gene-sRNA-target-GO networks
- **Platform**: Desktop GUI for Windows/macOS/Linux

---

## ✅ Summary Table

| Tool/Language     | Type          | Platform         | Purpose                                   |
|------------------|---------------|------------------|-------------------------------------------|
| Python 3.11       | Language      | Scripting/Modeling| Core analysis pipeline                    |
| R 4.4.1           | Language      | Scripting         | Enrichment, plots                         |
| miRDeep2          | Software      | Linux/CLI         | miRNA annotation                          |
| Bowtie            | Software      | Linux/CLI         | Read mapping                              |
| SPORTS1.1         | Software      | CLI               | tsRNA/rsRNA annotation                    |
| miRanda           | Software      | CLI               | Target gene prediction                    |
| RNAhybrid         | Software      | CLI               | Target binding energy calculation         |
| Cytoscape         | Desktop GUI   | Cross-platform    | Network visualization                     |

