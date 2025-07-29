# MDD-sRNA Diagnostic Modeling Framework

This project implements the full pipeline described in the paper "Circulating Small RNA signature for Diagnosis of Major Depressive Disorder", including small RNA expression processing, statistical analysis, HAMD-17 prediction modeling via Elastic net regression, and multi-objective model prioritization (MOMP)-based classifier design.  

## 1. Project Structure
```
MDD-sRNA-Diagnostic/
├── data/                         # Sample data. Fully preprocessed data matrices can be accessed via the Supplementary materials.
├── notebooks/                    # Jupyter notebooks for analysis
├── scripts&notebooks/                      # Python scripts or Jupyter notebooks for analysis
│   ├── preprocessing/            # Data filtering, normalization, annotation
│   ├── modeling/                 # MOMP optimization, ML model training
│   ├── validation/               # Cross-validation, AUC computation
│   ├── feature_selection/        # SHAP, clustering, forward selection
│   └── rt_qpcr_modeling/         # RT-qPCR validation and final scoring
├── results/                      # Output files, figures, and tables
├── models/                       # Trained models
├── references/                   # References and citations
├── README.md                     # Project overview
└── requirements.txt              # Python dependencies
```
# How to run  

To run it on your own computer: 
1. Install Python以及相关依赖项
3. Clone this repository
4. run the Jupyter Notbook Back_averaging_myoclonus.ipynb
  - from terminal: https://jupyter.readthedocs.io/en/latest/running.html
  - alternative for Mac users: https://nteract.io/desktop

# Dependencies
mne:           0.14
numpy:         1.11.3
scipy:         0.19.0
matplotlib:    1.5.1

sklearn:       0.18.1
nibabel:       2.1.0
mayavi:        4.5.0
pandas:        0.19.2
