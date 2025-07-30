# merge.py

import os
import pandas as pd

folder_path = './'
file_paths = [os.path.join(folder_path, f) for f in os.listdir(folder_path) if f.endswith('header.csv')]

merged_df = None

for file_path in file_paths:
    df = pd.read_csv(file_path, header=None, sep='\t')

    sample_name = os.path.basename(file_path).replace('_header.csv', '')

    if merged_df is None:
        df.columns = ['mirna_name_precursor', sample_name]
        merged_df = df
    else:
        df.columns = ['mirna_name_precursor', sample_name]
        merged_df = pd.merge(merged_df, df, on='mirna_name_precursor', how='outer')

merged_df = merged_df.fillna(0)
merged_df.to_csv('merged_mirna.csv', index=False, sep=',')
print("✅ Merged matrix saved to merged_mirna.csv")
