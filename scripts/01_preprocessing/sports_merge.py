import pandas as pd
import glob
import os

folder_path = './sports_merge'
file_paths = glob.glob(os.path.join(folder_path, '3*.csv'))

all_unique_subclass = set()

for file_path in file_paths:
    df = pd.read_csv(file_path)
    if 'Sequence:Annotation' in df.columns:
        all_unique_subclass.update(df['Sequence:Annotation'].tolist())

result_df = pd.DataFrame({'Sequence:Annotation': sorted(list(all_unique_subclass))})

for file_path in file_paths:
    identifier = file_path.split('_')[1]
    df = pd.read_csv(file_path)
    df = df.rename(columns={'Reads': f'Reads_{identifier}'})
    result_df = pd.merge(result_df, df, on='Sequence:Annotation', how='outer')

result_df = result_df.fillna(0)
result_df.to_csv('./sports_merge/merged_result.csv', index=False)
