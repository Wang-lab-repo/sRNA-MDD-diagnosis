# sRNA-MDD-diagnosis

Repository with scripts for all model training and analysis for paper "Circulating Small RNA signature for Diagnosis of Major Depressive Disorder"

The fully pre-processed input data for training the models can be found in the supplementary files. For each cancer, the basic data we used is **'CANCER_DATA_TOP2_JOINED_BATCH_CORRECTED_CLEANED.tsv'** where CANCER is the name of the cancer type. This data is GEO datasets collected from top 2 platforms, intersecting genes taken, and batch correction applied.

# Myoclonus back averaging

To run it on your own computer: 
1. Install Anaconda: https://docs.anaconda.com/anaconda/install/
2. Install MNE-Python: http://martinos.org/mne/stable/install_mne_python.html 
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
