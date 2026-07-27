"""
integrate_batches.py
====================
Batch integration and preprocessing of sRNA expression data.

Loads miRNA, rsRNA, and tsRNA expression matrices from two sequencing batches,
merges them with clinical metadata, applies CPM normalization, filters features
by expression prevalence (>=30% in ZMD+GZBrain), and splits into three cohorts:

  - ZMD_GZ.csv        : ZMD + GZBrain (training/internal validation)
  - NJBold_LYG.csv    : NJB_old + LYG (external validation 1)
  - NJBnew_Tower.csv  : NJB_new + Tower (external validation 2)

Expected input data structure:
  ./batch1/  : merged_mirna.csv, rsRNA_data.txt, tsRNA_data.txt, batch1_info.csv
  ./batch2/  : merged_mirna.csv, rsRNA_result.csv, tsRNA_result.csv, batch2_info.csv

Output:
  ./data/    : ZMD_GZ.csv, NJBold_LYG.csv, NJBnew_Tower.csv
"""
import pandas as pd
import numpy as np
import os
import warnings
warnings.filterwarnings('ignore')

# ================== Paths ==================
BATCH1_MIRNA = "./batch1/merged_mirna.csv"
BATCH1_RSRNA = "./batch1/rsRNA_data.txt"
BATCH1_TSRNA = "./batch1/tsRNA_data.txt"
BATCH1_CLINICAL = "./batch1/batch1_info.csv"

BATCH2_MIRNA = "./batch2/merged_mirna.csv"
BATCH2_RSRNA = "./batch2/rsRNA_result.csv"
BATCH2_TSRNA = "./batch2/tsRNA_result.csv"
BATCH2_CLINICAL = "./batch2/batch2_info.csv"

OUTPUT_DIR = "./data"
OUTPUT_ZMD_GZ = os.path.join(OUTPUT_DIR, "ZMD_GZ.csv")
OUTPUT_NJBOLD_LYG = os.path.join(OUTPUT_DIR, "NJBold_LYG.csv")
OUTPUT_NJBNEW_TOWER = os.path.join(OUTPUT_DIR, "NJBnew_Tower.csv")
MIN_SAMPLE_FRACTION = 0.3

# ================== Utility functions ==================
def read_expr_direct(filepath, prefix, sep='\t'):
    df = pd.read_csv(filepath, sep=sep, index_col=0)
    df = df.apply(pd.to_numeric, errors='coerce').fillna(0).astype(int)
    df.columns = [f"{prefix}_{c}" for c in df.columns]
    df.index.name = 'sample_id'
    return df

def cpm_normalize(expr_df):
    lib_size = expr_df.sum(axis=1)
    lib_size[lib_size == 0] = 1
    cpm = expr_df.div(lib_size, axis=0) * 1e6
    return cpm

def merge_rna_batch(mirna_df, rsrna_df, tsrna_df):
    for d in [mirna_df, rsrna_df, tsrna_df]:
        d.index = d.index.astype(str)
    merged = pd.concat([mirna_df, rsrna_df, tsrna_df], axis=1, join='outer')
    return merged.fillna(0)

# ================== 1. Load batch-1 ==================
print("Loading batch 1 ...")
mirna1 = pd.read_csv(BATCH1_MIRNA, index_col=0)
mirna1 = mirna1.apply(pd.to_numeric, errors='coerce').fillna(0).astype(int).T
mirna1.columns = [f"miRNA_{c}" for c in mirna1.columns]
mirna1.index.name = 'sample_id'

rsrna1 = read_expr_direct(BATCH1_RSRNA, 'rsRNA')
tsrna1 = read_expr_direct(BATCH1_TSRNA, 'tsRNA')

expr1 = cpm_normalize(merge_rna_batch(mirna1, rsrna1, tsrna1))

clinical1 = pd.read_csv(BATCH1_CLINICAL)
clinical1['Sample_id'] = clinical1.iloc[:, 1].astype(str)
clinical1['full_id'] = clinical1.iloc[:, 0].astype(str)
expr1.index = expr1.index.str.replace('@', '.')
expr1.index.name = 'full_id'
merged1 = expr1.merge(clinical1[['full_id', 'Sample_id']], left_index=True, right_on='full_id', how='inner')
merged1.set_index('Sample_id', inplace=True)
merged1.drop(columns=['full_id'], inplace=True, errors='ignore')
feature_cols1 = [c for c in merged1.columns if c.startswith(('miRNA_', 'rsRNA_', 'tsRNA_'))]
expr1_data = merged1[feature_cols1].astype(float)
obs1 = merged1.drop(columns=feature_cols1)
obs1['source_batch'] = 'batch1'

# ================== 2. Load batch-2 ==================
print("Loading batch 2 ...")
mirna2 = pd.read_csv(BATCH2_MIRNA, skiprows=[0], index_col=0)
mirna2 = mirna2.apply(pd.to_numeric, errors='coerce').fillna(0).astype(int).T
mirna2.columns = [f"miRNA_{c}" for c in mirna2.columns]
mirna2.index.name = 'sample_id'

rsrna2 = pd.read_csv(BATCH2_RSRNA, index_col=0)
rsrna2.index = rsrna2.index.str.split(':').str[0]
rsrna2 = rsrna2.apply(pd.to_numeric, errors='coerce').fillna(0).astype(int).T
rsrna2.columns = [f"rsRNA_{c}" for c in rsrna2.columns]

tsrna2 = pd.read_csv(BATCH2_TSRNA, index_col=0)
tsrna2.index = tsrna2.index.str.split(':').str[0]
tsrna2 = tsrna2.apply(pd.to_numeric, errors='coerce').fillna(0).astype(int).T
tsrna2.columns = [f"tsRNA_{c}" for c in tsrna2.columns]

expr2 = cpm_normalize(merge_rna_batch(mirna2, rsrna2, tsrna2))
common_rna2 = list(set(mirna2.index) & set(rsrna2.index) & set(tsrna2.index))
expr2 = expr2.loc[common_rna2]

clinical2 = pd.read_csv(BATCH2_CLINICAL)
clinical2['Sample_id'] = clinical2['Sample_id'].astype(str)
common_samp2 = list(set(expr2.index) & set(clinical2['Sample_id']))
expr2 = expr2.loc[common_samp2]
obs2 = clinical2.set_index('Sample_id').loc[common_samp2]
feature_cols2 = [c for c in expr2.columns if c.startswith(('miRNA_','rsRNA_','tsRNA_'))]
expr2_data = expr2[feature_cols2].astype(float)
obs2['source_batch'] = 'batch2'

# ================== 3. Merge batches (union of features, missing -> 0) ==================
print("Merging batches (union of features, missing -> 0) ...")
all_obs = pd.concat([obs1, obs2], axis=0, join='outer').fillna('NA')
all_expr = pd.concat([expr1_data, expr2_data], axis=0, join='outer').fillna(0)

# ================== 4. Assign NJB_old / NJB_new by original batch ==================
njb_mask = all_obs['Hospital'] == 'NJBrain'
all_obs.loc[njb_mask & (all_obs['source_batch'] == 'batch1'), 'Hospital'] = 'NJB_old'
all_obs.loc[njb_mask & (all_obs['source_batch'] == 'batch2'), 'Hospital'] = 'NJB_new'
print("NJB_old:", (all_obs['Hospital'] == 'NJB_old').sum(),
      "NJB_new:", (all_obs['Hospital'] == 'NJB_new').sum())

# ================== 5. Feature filtering on ZMD + GZBrain (>=30% expression) ==================
print("Filtering features using ZMD + GZBrain (>=30% samples with expression > 0) ...")
filter_mask = all_obs['Hospital'].isin(['ZMD', 'GZBrain'])
filter_expr = all_expr.loc[filter_mask]
prevalence = (filter_expr > 0).mean(axis=0)
filtered_features = prevalence[prevalence >= MIN_SAMPLE_FRACTION].index.tolist()
print(f"Features retained: {len(filtered_features)} / {all_expr.shape[1]}")

# ================== 6. Align all samples to filtered features ==================
all_expr_filtered = all_expr[filtered_features].fillna(0)

# ================== 7. Split and save ==================
os.makedirs(OUTPUT_DIR, exist_ok=True)

# ZMD + GZBrain
mask_zmd_gz = all_obs['Hospital'].isin(['ZMD', 'GZBrain'])
pd.concat([all_obs.loc[mask_zmd_gz], all_expr_filtered.loc[mask_zmd_gz]], axis=1).to_csv(OUTPUT_ZMD_GZ)
print(f"Saved {OUTPUT_ZMD_GZ}: {mask_zmd_gz.sum()} samples, {len(filtered_features)} features")

# NJB_old + LYG
mask_njbold_lyg = all_obs['Hospital'].isin(['NJB_old', 'LYG'])
pd.concat([all_obs.loc[mask_njbold_lyg], all_expr_filtered.loc[mask_njbold_lyg]], axis=1).to_csv(OUTPUT_NJBOLD_LYG)
print(f"Saved {OUTPUT_NJBOLD_LYG}: {mask_njbold_lyg.sum()} samples, {len(filtered_features)} features")

# NJB_new + Tower
mask_njbnew_tower = all_obs['Hospital'].isin(['NJB_new', 'Tower'])
pd.concat([all_obs.loc[mask_njbnew_tower], all_expr_filtered.loc[mask_njbnew_tower]], axis=1).to_csv(OUTPUT_NJBNEW_TOWER)
print(f"Saved {OUTPUT_NJBNEW_TOWER}: {mask_njbnew_tower.sum()} samples, {len(filtered_features)} features")

print("Done.")