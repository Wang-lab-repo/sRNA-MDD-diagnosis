import sys
import pandas as pd


if len(sys.argv) != 3:
    print("Usage: python median.py <input_file_path> <output_file_path>")
    sys.exit(1)

input_file_path = sys.argv[1]
output_file_path = sys.argv[2]

try:
    data = pd.read_csv(input_file_path)
except FileNotFoundError:
    print(f"Error: The file '{input_file_path}' was not found.")
    sys.exit(1)
except pd.errors.EmptyDataError:
    print(f"Error: The file '{input_file_path}' is empty.")
    sys.exit(1)
except pd.errors.ParserError:
    print(f"Error: There was a problem parsing '{input_file_path}'. Ensure it is a valid CSV file.")
    sys.exit(1)


data = pd.read_csv(input_file_path)
mirna_names = data['miRNA']
control_columns = data.filter(regex='^control_')
treatment_columns = data.filter(regex='^treatment_')
control_median = control_columns.replace(0, pd.NA).median(axis=1, skipna=True)
treatment_median = treatment_columns.replace(0, pd.NA).median(axis=1, skipna=True)
fold_change = treatment_median - control_median
result_data = pd.DataFrame({'miRNA': mirna_names, 'log2control_median': control_median, 'log2treatment_median': treatment_median, 'log2fold_change': fold_change})
try:
    result_data.to_csv(output_file_path, index=False)
except Exception as e:
    print(f"Error: Could not write to file '{output_file_path}'. Details: {e}")
    sys.exit(1)
print(f"Results successfully saved to '{output_file_path}'")
