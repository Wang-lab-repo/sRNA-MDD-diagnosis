import sys
import pandas as pd

# 检查命令行参数是否正确
if len(sys.argv) != 3:
    print("Usage: python median.py <input_file_path> <output_file_path>")
    sys.exit(1)

# 从命令行参数获取输入文件路径和输出文件路径
input_file_path = sys.argv[1]
output_file_path = sys.argv[2]

try:
    # 读取数据
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


# 读取数据
data = pd.read_csv(input_file_path)
# 提取 miRNA 名称列和相关的 control 和 treatment 列
mirna_names = data['rsRNA']
control_columns = data.filter(regex='^control_')
treatment_columns = data.filter(regex='^treatment_')
# 计算中位数
control_median = control_columns.replace(0, pd.NA).median(axis=1, skipna=True)
treatment_median = treatment_columns.replace(0, pd.NA).median(axis=1, skipna=True)
# 计算 Fold Change
fold_change = treatment_median - control_median
# 创建包含结果的新 DataFrame
result_data = pd.DataFrame({'rsRNA': mirna_names, 'log2control_median': control_median, 'log2treatment_median': treatment_median, 'log2fold_change': fold_change})
# 将结果保存到 CSV 文件
try:
    result_data.to_csv(output_file_path, index=False)
except Exception as e:
    print(f"Error: Could not write to file '{output_file_path}'. Details: {e}")
    sys.exit(1)
print(f"Results successfully saved to '{output_file_path}'")