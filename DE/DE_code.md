# 2. pirna差异分析  
  
1. 统计每个样本中表达量不为0的rna种类
```r
library("ggplot2")
setwd("D:/adult_dep/DE2/")
df <- read.csv("./input/rpm_50%_log_pirna.csv", header = TRUE, row.names = 1)
df <- df[, df[1, ] %in% c(0, 1)]
df <- df[-1, ]
result_df <- data.frame()

for (sample_col in 2:ncol(df)) {
  sample_name <- colnames(df)[sample_col]
  non_zero_mirna_count <- sum(df[, sample_col] != 0)
  result_df <- rbind(result_df, data.frame(sample = sample_name, non_zero_mirna_count = non_zero_mirna_count))
}
print(result_df)
write.csv(result_df, file = "./MDD_hc_pirna/non_zero_pirna_count.csv")
mean_non_zero_mirna_count <- mean(result_df$non_zero_mirna_count)
cat("Mean of non_zero_mirna_count:", mean_non_zero_mirna_count, "\n")
# 平均样本的mirna数目为37.66507 

ggplot(result_df, aes(x = non_zero_mirna_count)) +
  geom_histogram(binwidth = 1, fill = "blue", color = "black", alpha = 0.7) +
  labs(title = "Frequency Distribution of non_zero_mirna_count",
       x = "non_zero_mirna_count",
       y = "Frequency")
```

2. 标注control_vs_treatment，手动标注
```r
df <- read.csv("./input/rpm_50%_log_pirna.csv", header = TRUE, row.names = 1)
# 只保留第一行中 group 为 0 或 1 的列
df <- df[, df[1, ] %in% c(0, 1)]
# 根据第一行的值（0 或 1）修改列名
colnames(df) <- ifelse(df[1, ] == 0, 
                       paste0("control_", colnames(df)), 
                       paste0("treatment_", colnames(df)))
# 删除第一行
df <- df[-1, ]
# 查看修改后的数据框
head(df)
write.csv(df, file = "./input_rename/pirna.csv")
```

3. 找到每个pirna的treatment_vs_control中位数，计算FC
```bash
cd /mnt/d/adult_dep/DE2
# 记得加“pirna”
python median.py ./input_rename/pirna.csv ./MDD_hc_pirna/median.csv
```
```python
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
mirna_names = data['pirna']
control_columns = data.filter(regex='^control_')
treatment_columns = data.filter(regex='^treatment_')
# 计算中位数
control_median = control_columns.replace(0, pd.NA).median(axis=1, skipna=True)
treatment_median = treatment_columns.replace(0, pd.NA).median(axis=1, skipna=True)
# 计算 Fold Change
fold_change = treatment_median - control_median
# 创建包含结果的新 DataFrame
result_data = pd.DataFrame({'piRNA': mirna_names, 'log2control_median': control_median, 'log2treatment_median': treatment_median, 'log2fold_change': fold_change})
# 将结果保存到 CSV 文件
try:
    result_data.to_csv(output_file_path, index=False)
except Exception as e:
    print(f"Error: Could not write to file '{output_file_path}'. Details: {e}")
    sys.exit(1)
print(f"Results successfully saved to '{output_file_path}'")
```


② 缺失值替换为NA
```r
setwd("D:/adult_dep/DE2")
data <- read.csv("./input_rename/pirna.csv") 
data[data == 0] <- NA
write.csv(data, "./input_rename/pirna_NA.csv", row.names = FALSE)  
```



③ Wilcoxon秩和检验（相关样本）  

秩和检验（Wilcoxon秩和检验）是一种用于比较两组相关或无关样本的非参数检验方法。它不要求数据符合正态分布，因此适用于那些无法满足正态性假设的情况。Wilcoxon秩和检验（相关样本）： 适用于配对设计，比如同一组个体在两个不同时间点的观测。在这种情况下，对每对相关样本进行秩和的比较。    
* 原理：在Wilcoxon秩和检验中，基本思想是将每个组的观测值按大小排序，然后为每个观测值分配一个秩次，最后通过比较秩和来判断两组样本是否来自同一分布。这个检验的零假设是两组样本来自相同的分布。  


```r
# 设置工作目录
setwd("D:/adult_dep/DE2")

# 读取数据
data <- read.csv("./input_rename/pirna_NA.csv")

# 加载dplyr包
library(dplyr)

# 获取控制组和治疗组的列
control_columns <- grep("^control_", names(data), value = TRUE)
treatment_columns <- grep("^treatment_", names(data), value = TRUE)
# 835

# 执行 Wilcoxon 秩和检验
p_values_matrix <- matrix(NA, nrow = length(data$piRNA), ncol = 3)
colnames(p_values_matrix) <- c("piRNA", "p_value", "Q_value")

for (i in seq_along(data$piRNA)) {
  piRNA <- data$piRNA[i]
  control_values <- unlist(data[data$piRNA == piRNA, control_columns])
  treatment_values <- unlist(data[data$piRNA == piRNA, treatment_columns])
  
  # 执行 Wilcoxon 秩和检验
  result <- try(wilcox.test(treatment_values, control_values))

  # 存储 p-value
  if (!inherits(result, "try-error")) {
    p_values_matrix[i, 1] <- piRNA
    p_values_matrix[i, 2] <- result$p.value
  } else {
    cat("Error for piRNA:", piRNA, "\n")
  }
}

# 使用 Benjamini-Hochberg 方法进行多重比较校正
p_values_matrix[, 3] <- p.adjust(p_values_matrix[, 2], method = "BH")

# 输出结果
print(p_values_matrix)
# 将 p_values_matrix 保存到 CSV 文件
write.csv(p_values_matrix, file = "./MDD_hc_pirna/de_values_matrix.csv", row.names = FALSE)
```

④ 合并文件
```r
library(dplyr)

median_data <- read.csv("./MDD_hc_pirna/median.csv")
de_values_data <- read.csv("./MDD_hc_pirna/de_values_matrix.csv")
# 基于 miRNA 列合并数据
merged_data <- inner_join(median_data, de_values_data, by = "piRNA")

# 查看合并后的数据
print(merged_data)
write.csv(merged_data, file = "./MDD_hc_pirna/result.csv", row.names = FALSE)

# 4 up, 19 down

# padj 小于 0.05 并且 Log2FC 大于 1（2倍） 或者小于 -1（1/2倍）
diff_gene <- subset(merged_data, Q_value < 0.05 & abs(log2fold_change) > 1)
write.csv(diff_gene, file="./MDD_hc_pirna/sig_q0.05_fc2.csv", quote = F)
```
```r
library(dplyr)
setwd("D:/adult_dep/DE2")

merged_data <- read.csv("./MDD_hc_pirna/result.csv")

diff_gene <- subset(merged_data, Q_value < 0.05 & abs(log2fold_change) > 0.58)
write.csv(diff_gene, file="./MDD_hc_pirna/sig_q0.05_fc15.csv", quote = F)
```


⑤ 绘图
```r
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}
if (!requireNamespace("pheatmap", quietly = TRUE)) {
  install.packages("pheatmap")
}

# 加载包
library(ggplot2)
library(pheatmap)

merged_data <- read.csv("./MDD_hc_pirna/result.csv")
merged_data$lgpadj <- -log10(merged_data$Q_value)

# 绘制火山图
volcano_plot <- ggplot(merged_data, aes(x = log2fold_change, y = lgpadj)) +
  geom_point(aes(color = ifelse(Q_value < 0.05 & log2fold_change > 0.58, "Up", 
                                ifelse(Q_value < 0.05 & log2fold_change < -0.58, "Down", "Not Sig"))), size = 2) +
  scale_color_manual(values = c("Up" = "#D32F2F", "Down" = "#1976D2", "Not Sig" = "grey70")) +
  labs(x = "Log2 Fold Change", y = "-Log10(padj)", title = "Volcano Plot") +
  theme_minimal(base_size = 14) +
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(),  # 去掉主网格线
        panel.grid.minor = element_blank(),  # 去掉次网格线
        axis.line = element_line(color = "black"), # 显示横纵坐标轴
        axis.ticks = element_line(color = "black"),  # 设置坐标轴刻度线颜色
        axis.ticks.length = unit(0.25, "cm")) +  # 设置坐标轴刻度线长度
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", size = 0.5) +
  geom_vline(xintercept = c(-0.58, 0.58), linetype = "dashed", color = "black", size = 0.5)

ggsave("./MDD_hc_pirna/MDD_hc_pirna_volcano_plot.pdf", plot = volcano_plot, device = "pdf", dpi = 600, width = 8, height = 6)
```
* 图中标注

```r
# 图中标注名称
miRNAs_to_annotate <- c("piR-hsa-143928", "piR-hsa-369497")

# 提取需要标注的 miRNA 数据
annotations <- merged_data[merged_data$piRNA %in% miRNAs_to_annotate, ]

# 创建火山图
volcano_plot <- ggplot(merged_data, aes(x = log2fold_change, y = lgpadj)) +
  geom_point(aes(color = ifelse(Q_value < 0.05 & log2fold_change > 0.58, "Up", 
                                ifelse(Q_value < 0.05 & log2fold_change < -0.58, "Down", "Not Sig"))), size = 2) +
  scale_color_manual(values = c("Up" = "#D32F2F", "Down" = "#1976D2", "Not Sig" = "grey70")) +
  labs(x = "Log2 Fold Change", y = "-Log10(padj)", title = "Volcano Plot") +
  theme_minimal(base_size = 14) +
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", size = 0.5) +
  geom_vline(xintercept = c(-0.58, 0.58), linetype = "dashed", color = "black", size = 0.5) +
  geom_point(data = annotations, aes(color = "Highlighted"), size = 3, shape = 1, stroke = 1.5, fill = NA) +  # 添加圈圈
  geom_text(data = annotations, aes(label = piRNA), vjust = -1, hjust = 0.5, size = 3, color = "black")


ggsave("./MDD_hc_pirna/label_MDD_hc_pirna_volcano_plot.pdf", plot = volcano_plot, device = "pdf", dpi = 600, width = 8, height = 6)

```

```r
# 绘制热图
library(ggplot2)
library(tidyr)
library(dplyr)

# 读取数据
merged_data <- read.csv("./MDD_hc_pirna/result.csv")

# 指定要绘制热图的 miRNA
piRNAs_to_annotate <- c("piR-hsa-143928", "piR-hsa-369497")

# 提取需要的 piRNA 数据并设置顺序
heatmap_data <- merged_data %>%
  filter(piRNA %in% piRNAs_to_annotate) %>%
  select(piRNA, log2control_median, log2treatment_median) %>%
  pivot_longer(cols = starts_with("log2"), names_to = "Condition", values_to = "Expression") %>%
  mutate(piRNA = factor(piRNA, levels = piRNAs_to_annotate))  # 设置因子顺序

# 定义颜色从白色到深蓝色
start_color <- rgb(255/255, 255/255, 255/255)  # 白色
end_color <- rgb(0/255, 51/255, 102/255)        # 深蓝色

# 生成渐变
num_colors <- 10
gradient_colors <- colorRampPalette(c(start_color, end_color))(num_colors)

# 绘制热图
heatmap_plot <- ggplot(heatmap_data, aes(x = piRNA, y = Condition, fill = Expression)) +
  geom_tile(color = "white") +
  scale_fill_gradientn(colors = gradient_colors, 
                       name = "Expression Level") +
  labs(x = "piRNA", y = "Condition", title = "Heatmap of piRNA Expression") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 保存热图
ggsave("./MDD_hc_piRNA/5piRNA_heatmap.pdf", plot = heatmap_plot, device = "pdf", dpi = 600, width = 8, height = 6)
```


















# 3. tsrna差异分析  
  
1. 统计每个样本中表达量不为0的rna种类
```r
library("ggplot2")
setwd("D:/adult_dep/DE2/")
df <- read.csv("./input/rpm_50%_log_tsrna.csv", header = TRUE, row.names = 1)
df <- df[, df[1, ] %in% c(0, 1)]
df <- df[-1, ]
result_df <- data.frame()

for (sample_col in 2:ncol(df)) {
  sample_name <- colnames(df)[sample_col]
  non_zero_mirna_count <- sum(df[, sample_col] != 0)
  result_df <- rbind(result_df, data.frame(sample = sample_name, non_zero_mirna_count = non_zero_mirna_count))
}
print(result_df)
write.csv(result_df, file = "./MDD_hc_tsrna/non_zero_tsrna_count.csv")
mean_non_zero_mirna_count <- mean(result_df$non_zero_mirna_count)
cat("Mean of non_zero_tsrna_count:", mean_non_zero_mirna_count, "\n")
# 平均样本的tsrna数目为62.44656/99

ggplot(result_df, aes(x = non_zero_mirna_count)) +
  geom_histogram(binwidth = 1, fill = "blue", color = "black", alpha = 0.7) +
  labs(title = "Frequency Distribution of non_zero_mirna_count",
       x = "non_zero_mirna_count",
       y = "Frequency")
```

2. 标注control_vs_treatment，手动标注
```r
df <- read.csv("./input/rpm_50%_log_tsrna.csv", header = TRUE, row.names = 1)
# 只保留第一行中 group 为 0 或 1 的列
df <- df[, df[1, ] %in% c(0, 1)]
# 根据第一行的值（0 或 1）修改列名
colnames(df) <- ifelse(df[1, ] == 0, 
                       paste0("control_", colnames(df)), 
                       paste0("treatment_", colnames(df)))
# 删除第一行
df <- df[-1, ]
# 查看修改后的数据框
head(df)
write.csv(df, file = "./input_rename/tsrna.csv")
```

3. 找到每个pirna的treatment_vs_control中位数，计算FC
```bash
cd /mnt/d/adult_dep/DE2
# 记得加“tsrna”
python median.py ./input_rename/tsrna.csv ./MDD_hc_tsrna/median.csv
```
```python
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
mirna_names = data['tsRNA']
control_columns = data.filter(regex='^control_')
treatment_columns = data.filter(regex='^treatment_')
# 计算中位数
control_median = control_columns.replace(0, pd.NA).median(axis=1, skipna=True)
treatment_median = treatment_columns.replace(0, pd.NA).median(axis=1, skipna=True)
# 计算 Fold Change
fold_change = treatment_median - control_median
# 创建包含结果的新 DataFrame
result_data = pd.DataFrame({'tsRNA': mirna_names, 'log2control_median': control_median, 'log2treatment_median': treatment_median, 'log2fold_change': fold_change})
# 将结果保存到 CSV 文件
try:
    result_data.to_csv(output_file_path, index=False)
except Exception as e:
    print(f"Error: Could not write to file '{output_file_path}'. Details: {e}")
    sys.exit(1)
print(f"Results successfully saved to '{output_file_path}'")
```


② 缺失值替换为NA
```r
setwd("D:/adult_dep/DE2")
data <- read.csv("./input_rename/tsRNA.csv") 
data[data == 0] <- NA
write.csv(data, "./input_rename/tsRNA_NA.csv", row.names = FALSE)  
```



③ Wilcoxon秩和检验（相关样本）  

秩和检验（Wilcoxon秩和检验）是一种用于比较两组相关或无关样本的非参数检验方法。它不要求数据符合正态分布，因此适用于那些无法满足正态性假设的情况。Wilcoxon秩和检验（相关样本）： 适用于配对设计，比如同一组个体在两个不同时间点的观测。在这种情况下，对每对相关样本进行秩和的比较。    
* 原理：在Wilcoxon秩和检验中，基本思想是将每个组的观测值按大小排序，然后为每个观测值分配一个秩次，最后通过比较秩和来判断两组样本是否来自同一分布。这个检验的零假设是两组样本来自相同的分布。  


```r
# 设置工作目录
setwd("D:/adult_dep/DE2")

# 读取数据
data <- read.csv("./input_rename/tsrna_NA.csv")

# 加载dplyr包
library(dplyr)

# 获取控制组和治疗组的列
control_columns <- grep("^control_", names(data), value = TRUE)
treatment_columns <- grep("^treatment_", names(data), value = TRUE)
# 835

# 执行 Wilcoxon 秩和检验
p_values_matrix <- matrix(NA, nrow = length(data$tsRNA), ncol = 3)
colnames(p_values_matrix) <- c("tsRNA", "p_value", "Q_value")

for (i in seq_along(data$tsRNA)) {
  tsRNA <- data$tsRNA[i]
  control_values <- unlist(data[data$tsRNA == tsRNA, control_columns])
  treatment_values <- unlist(data[data$tsRNA == tsRNA, treatment_columns])
  
  # 执行 Wilcoxon 秩和检验
  result <- try(wilcox.test(treatment_values, control_values))

  # 存储 p-value
  if (!inherits(result, "try-error")) {
    p_values_matrix[i, 1] <- tsRNA
    p_values_matrix[i, 2] <- result$p.value
  } else {
    cat("Error for tsRNA:", tsRNA, "\n")
  }
}

# 使用 Benjamini-Hochberg 方法进行多重比较校正
p_values_matrix[, 3] <- p.adjust(p_values_matrix[, 2], method = "BH")

# 输出结果
print(p_values_matrix)
# 将 p_values_matrix 保存到 CSV 文件
write.csv(p_values_matrix, file = "./MDD_hc_tsrna/de_values_matrix.csv", row.names = FALSE)
```

④ 合并文件
```r
library(dplyr)

median_data <- read.csv("./MDD_hc_tsrna/median.csv")
de_values_data <- read.csv("./MDD_hc_tsrna/de_values_matrix.csv")
# 基于 miRNA 列合并数据
merged_data <- inner_join(median_data, de_values_data, by = "tsRNA")

# 查看合并后的数据
print(merged_data)
write.csv(merged_data, file = "./MDD_hc_tsrna/result.csv", row.names = FALSE)

# padj 小于 0.05 并且 Log2FC 大于 1（2倍） 或者小于 -1（1/2倍）
diff_gene <- subset(merged_data, Q_value < 0.05 & abs(log2fold_change) > 1)
write.csv(diff_gene, file="./MDD_hc_tsrna/sig_q0.05_fc2.csv", quote = F)
```

```r
library(dplyr)
setwd("D:/adult_dep/DE2")

merged_data <- read.csv("./MDD_hc_tsrna/result.csv")

diff_gene <- subset(merged_data, Q_value < 0.05 & abs(log2fold_change) > 0.58)
write.csv(diff_gene, file="./MDD_hc_tsrna/sig_q0.05_fc15.csv", quote = F)
```


⑤ 绘图
```r
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}
if (!requireNamespace("pheatmap", quietly = TRUE)) {
  install.packages("pheatmap")
}

# 加载包
library(ggplot2)
library(pheatmap)

merged_data <- read.csv("./MDD_hc_tsrna/result.csv")
merged_data$lgpadj <- -log10(merged_data$Q_value)

# 绘制火山图
volcano_plot <- ggplot(merged_data, aes(x = log2fold_change, y = lgpadj)) +
  geom_point(aes(color = ifelse(Q_value < 0.05 & log2fold_change > 0.58, "Up", 
                                ifelse(Q_value < 0.05 & log2fold_change < -0.58, "Down", "Not Sig"))), size = 2) +
  scale_color_manual(values = c("Up" = "#D32F2F", "Down" = "#1976D2", "Not Sig" = "grey70")) +
  labs(x = "Log2 Fold Change", y = "-Log10(padj)", title = "Volcano Plot") +
  theme_minimal(base_size = 14) +
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(),  # 去掉主网格线
        panel.grid.minor = element_blank(),  # 去掉次网格线
        axis.line = element_line(color = "black"), # 显示横纵坐标轴
        axis.ticks = element_line(color = "black"),  # 设置坐标轴刻度线颜色
        axis.ticks.length = unit(0.25, "cm")) +  # 设置坐标轴刻度线长度
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", size = 0.5) +
  geom_vline(xintercept = c(-0.58, 0.58), linetype = "dashed", color = "black", size = 0.5)

ggsave("./MDD_hc_tsrna/MDD_hc_tsrna_volcano_plot.pdf", plot = volcano_plot, device = "pdf", dpi = 600, width = 8, height = 6)
```
* 图中标注

```r
# 图中标注名称
annotations <- data.frame(
  tsRNA = c("TTGGTCGTGGTTGTAGTCCGTGCGAGAATACCA", "ATCCCGGACGAGCCCCCA", 
            "CCCTGGTGGTCTAGTGGTTAGGATTCGGCGCT", "GCATTGGTGGTTCAGTGGTAGAATTCTCGCCT"),
  new_label = c("tsRNA-#98", "tsRNA-#16", "tsRNA-#25", "tsRNA-#41")
)

# 提取需要标注的 miRNA 数据
annotations_data <- merged_data[merged_data$tsRNA %in% annotations$tsRNA, ]

# 创建火山图
volcano_plot <- ggplot(merged_data, aes(x = log2fold_change, y = lgpadj)) +
  geom_point(aes(color = ifelse(Q_value < 0.05 & log2fold_change > 0.58, "Up", 
                                ifelse(Q_value < 0.05 & log2fold_change < -0.58, "Down", "Not Sig"))), size = 2) +
  scale_color_manual(values = c("Up" = "#D32F2F", "Down" = "#1976D2", "Not Sig" = "grey70")) +
  labs(x = "Log2 Fold Change", y = "-Log10(padj)", title = "Volcano Plot") +
  theme_minimal(base_size = 14) +
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", size = 0.5) +
  geom_vline(xintercept = c(-0.58, 0.58), linetype = "dashed", color = "black", size = 0.5) +
  geom_point(data = annotations_data, aes(color = "Highlighted"), size = 3, shape = 1, stroke = 1.5, fill = NA) +  # 添加圈圈
  geom_text(data = merge(annotations_data, annotations, by = "tsRNA"), aes(label = new_label), 
            vjust = -1, hjust = 0.5, size = 3, color = "black")  # 使用新的标签


ggsave("./MDD_hc_tsrna/label_MDD_hc_tsrna_volcano_plot.pdf", plot = volcano_plot, device = "pdf", dpi = 600, width = 8, height = 6)

```

```r
# 绘制热图
library(ggplot2)
library(tidyr)
library(dplyr)

# 读取数据
merged_data <- read.csv("./MDD_hc_tsrna/result.csv")

# 指定需要标注的 tsRNA 和对应的新名称
tsRNAs_to_annotate <- c("TTGGTCGTGGTTGTAGTCCGTGCGAGAATACCA", "ATCCCGGACGAGCCCCCA", 
                         "CCCTGGTGGTCTAGTGGTTAGGATTCGGCGCT", "GCATTGGTGGTTCAGTGGTAGAATTCTCGCCT")

new_names <- c("tsRNA-#98", "tsRNA-#16", "tsRNA-#25", "tsRNA-#41")  # 新名称

# 创建一个命名映射
name_mapping <- setNames(new_names, tsRNAs_to_annotate)

# 提取需要的 tsRNA 数据并设置顺序
heatmap_data <- merged_data %>%
  filter(tsRNA %in% tsRNAs_to_annotate) %>%
  select(tsRNA, log2control_median, log2treatment_median) %>%
  pivot_longer(cols = starts_with("log2"), names_to = "Condition", values_to = "Expression") %>%
  mutate(tsRNA = factor(name_mapping[tsRNA], levels = new_names))  # 使用新名称

# 定义颜色从白色到深蓝色
start_color <- rgb(255/255, 255/255, 255/255)  # 白色
end_color <- rgb(0/255, 51/255, 102/255)        # 深蓝色

# 生成渐变
num_colors <- 10
gradient_colors <- colorRampPalette(c(start_color, end_color))(num_colors)

# 绘制热图
heatmap_plot <- ggplot(heatmap_data, aes(x = tsRNA, y = Condition, fill = Expression)) +
  geom_tile(color = "white") +
  scale_fill_gradientn(colors = gradient_colors, 
                       name = "Expression Level") +
  labs(x = "tsRNA", y = "Condition", title = "Heatmap of tsRNA Expression") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 保存热图
ggsave("./MDD_hc_tsRNA/5tsRNA_heatmap.pdf", plot = heatmap_plot, device = "pdf", dpi = 600, width = 8, height = 6)
```




























# 4. rsrna差异分析  
  
1. 统计每个样本中表达量不为0的rna种类
```r
library("ggplot2")
setwd("D:/adult_dep/DE2/")
df <- read.csv("./input/rpm_50%_log_rsrna.csv", header = TRUE, row.names = 1)
df <- df[, df[1, ] %in% c(0, 1)]
df <- df[-1, ]
result_df <- data.frame()

for (sample_col in 2:ncol(df)) {
  sample_name <- colnames(df)[sample_col]
  non_zero_mirna_count <- sum(df[, sample_col] != 0)
  result_df <- rbind(result_df, data.frame(sample = sample_name, non_zero_mirna_count = non_zero_mirna_count))
}
print(result_df)
write.csv(result_df, file = "./MDD_hc_rsrna/non_zero_rsrna_count.csv")
mean_non_zero_mirna_count <- mean(result_df$non_zero_mirna_count)
cat("Mean of non_zero_mirna_count:", mean_non_zero_mirna_count, "\n")
# 平均样本的mirna数目为37.66507 

ggplot(result_df, aes(x = non_zero_mirna_count)) +
  geom_histogram(binwidth = 1, fill = "blue", color = "black", alpha = 0.7) +
  labs(title = "Frequency Distribution of non_zero_mirna_count",
       x = "non_zero_mirna_count",
       y = "Frequency")
```

2. 标注control_vs_treatment，手动标注
```r
df <- read.csv("./input/rpm_50%_log_rsrna.csv", header = TRUE, row.names = 1)
# 只保留第一行中 group 为 0 或 1 的列
df <- df[, df[1, ] %in% c(0, 1)]
# 根据第一行的值（0 或 1）修改列名
colnames(df) <- ifelse(df[1, ] == 0, 
                       paste0("control_", colnames(df)), 
                       paste0("treatment_", colnames(df)))
# 删除第一行
df <- df[-1, ]
# 查看修改后的数据框
head(df)
write.csv(df, file = "./input_rename/rsrna.csv")
```

3. 找到每个rsrna的treatment_vs_control中位数，计算FC
```bash
cd /mnt/d/adult_dep/DE2
# 记得加“rsRNA”
python median.py ./input_rename/rsrna.csv ./MDD_hc_rsrna/median.csv
```
```python
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
```


② 缺失值替换为NA
```r
setwd("D:/adult_dep/DE2")
data <- read.csv("./input_rename/rsrna.csv") 
data[data == 0] <- NA
write.csv(data, "./input_rename/rsrna_NA.csv", row.names = FALSE)  
```



③ Wilcoxon秩和检验（相关样本）  

秩和检验（Wilcoxon秩和检验）是一种用于比较两组相关或无关样本的非参数检验方法。它不要求数据符合正态分布，因此适用于那些无法满足正态性假设的情况。Wilcoxon秩和检验（相关样本）： 适用于配对设计，比如同一组个体在两个不同时间点的观测。在这种情况下，对每对相关样本进行秩和的比较。    
* 原理：在Wilcoxon秩和检验中，基本思想是将每个组的观测值按大小排序，然后为每个观测值分配一个秩次，最后通过比较秩和来判断两组样本是否来自同一分布。这个检验的零假设是两组样本来自相同的分布。  


```r
# 设置工作目录
setwd("D:/adult_dep/DE2")

# 读取数据
data <- read.csv("./input_rename/rsrna_NA.csv")

# 加载dplyr包
library(dplyr)

# 获取控制组和治疗组的列
control_columns <- grep("^control_", names(data), value = TRUE)
treatment_columns <- grep("^treatment_", names(data), value = TRUE)
# 835

# 执行 Wilcoxon 秩和检验
p_values_matrix <- matrix(NA, nrow = length(data$rsRNA), ncol = 3)
colnames(p_values_matrix) <- c("rsRNA", "p_value", "Q_value")

for (i in seq_along(data$rsRNA)) {
  rsRNA <- data$rsRNA[i]
  control_values <- unlist(data[data$rsRNA == rsRNA, control_columns])
  treatment_values <- unlist(data[data$rsRNA == rsRNA, treatment_columns])
  
  # 执行 Wilcoxon 秩和检验
  result <- try(wilcox.test(treatment_values, control_values))

  # 存储 p-value
  if (!inherits(result, "try-error")) {
    p_values_matrix[i, 1] <- rsRNA
    p_values_matrix[i, 2] <- result$p.value
  } else {
    cat("Error for rsRNA:", rsRNA, "\n")
  }
}

# 使用 Benjamini-Hochberg 方法进行多重比较校正
p_values_matrix[, 3] <- p.adjust(p_values_matrix[, 2], method = "BH")

# 输出结果
print(p_values_matrix)
# 将 p_values_matrix 保存到 CSV 文件
write.csv(p_values_matrix, file = "./MDD_hc_rsrna/de_values_matrix.csv", row.names = FALSE)
```

④ 合并文件
```r
library(dplyr)

median_data <- read.csv("./MDD_hc_rsrna/median.csv")
de_values_data <- read.csv("./MDD_hc_rsrna/de_values_matrix.csv")
# 基于 miRNA 列合并数据
merged_data <- inner_join(median_data, de_values_data, by = "rsRNA")

# 查看合并后的数据
print(merged_data)
write.csv(merged_data, file = "./MDD_hc_rsrna/result.csv", row.names = FALSE)

# padj 小于 0.05 并且 Log2FC 大于 1（2倍） 或者小于 -1（1/2倍）
diff_gene <- subset(merged_data, Q_value < 0.05 & abs(log2fold_change) > 1)
write.csv(diff_gene, file="./MDD_hc_rsrna/sig_q0.05_fc2.csv", quote = F)
```

```r
library(dplyr)
setwd("D:/adult_dep/DE2")

merged_data <- read.csv("./MDD_hc_rsrna/result.csv")

diff_gene <- subset(merged_data, Q_value < 0.05 & abs(log2fold_change) > 0.58)
write.csv(diff_gene, file="./MDD_hc_rsrna/sig_q0.05_fc15.csv", quote = F)
```


⑤ 绘图
```r
# 加载包
library(ggplot2)
library(pheatmap)

merged_data <- read.csv("./MDD_hc_rsrna/result.csv")
merged_data$lgpadj <- -log10(merged_data$Q_value)

# 绘制火山图
volcano_plot <- ggplot(merged_data, aes(x = log2fold_change, y = lgpadj)) +
  geom_point(aes(color = ifelse(Q_value < 0.05 & log2fold_change > 0.58, "Up", 
                                ifelse(Q_value < 0.05 & log2fold_change < -0.58, "Down", "Not Sig"))), size = 2) +
  scale_color_manual(values = c("Up" = "#D32F2F", "Down" = "#1976D2", "Not Sig" = "grey70")) +
  labs(x = "Log2 Fold Change", y = "-Log10(padj)", title = "Volcano Plot") +
  theme_minimal(base_size = 14) +
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(),  # 去掉主网格线
        panel.grid.minor = element_blank(),  # 去掉次网格线
        axis.line = element_line(color = "black"), # 显示横纵坐标轴
        axis.ticks = element_line(color = "black"),  # 设置坐标轴刻度线颜色
        axis.ticks.length = unit(0.25, "cm")) +  # 设置坐标轴刻度线长度
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", size = 0.5) +
  geom_vline(xintercept = c(-0.58, 0.58), linetype = "dashed", color = "black", size = 0.5)

ggsave("./MDD_hc_rsrna/MDD_hc_rsrna_volcano_plot.pdf", plot = volcano_plot, device = "pdf", dpi = 600, width = 8, height = 6)
```
* 图中标注

```r
# 图中标注名称
annotations <- data.frame(
  rsRNA = c("CGGACCAAGGAGTCTAACACGTGCGCG", "GACCAAGGAGTCTAACACGTGCGCG", 
            "GTGGTGCATGGCCGTTCTT", "GTGGTGGTGCATGGCCGTTCTT", "TCGGTCGGGCTGGGGCGCGA"),
  new_label = c("rsRNA-#122", "rsRNA-#146", "rsRNA-#216", "rsRNA-#217", "rsRNA-#250")
)

# 提取需要标注的 miRNA 数据
annotations_data <- merged_data[merged_data$rsRNA %in% annotations$rsRNA, ]

# 创建火山图
volcano_plot <- ggplot(merged_data, aes(x = log2fold_change, y = lgpadj)) +
  geom_point(aes(color = ifelse(Q_value < 0.05 & log2fold_change > 0.58, "Up", 
                                ifelse(Q_value < 0.05 & log2fold_change < -0.58, "Down", "Not Sig"))), size = 2) +
  scale_color_manual(values = c("Up" = "#D32F2F", "Down" = "#1976D2", "Not Sig" = "grey70")) +
  labs(x = "Log2 Fold Change", y = "-Log10(padj)", title = "Volcano Plot") +
  theme_minimal(base_size = 14) +
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", size = 0.5) +
  geom_vline(xintercept = c(-0.58, 0.58), linetype = "dashed", color = "black", size = 0.5) +
  geom_point(data = annotations_data, aes(color = "Highlighted"), size = 3, shape = 1, stroke = 1.5, fill = NA) +  # 添加圈圈
  geom_text(data = merge(annotations_data, annotations, by = "rsRNA"), aes(label = new_label), 
            vjust = -1, hjust = 0.5, size = 3, color = "black")  # 使用新的标签


ggsave("./MDD_hc_rsrna/label_MDD_hc_rsrna_volcano_plot.pdf", plot = volcano_plot, device = "pdf", dpi = 600, width = 8, height = 6)

```

```r
# 绘制热图
library(ggplot2)
library(tidyr)
library(dplyr)
library(ggplot2)
library(pheatmap)

# 读取数据
merged_data <- read.csv("./MDD_hc_rsrna/result.csv")

# 指定需要标注的 rsRNA 和对应的新名称
rsRNAs_to_annotate <- c("CGGACCAAGGAGTCTAACACGTGCGCG", "GACCAAGGAGTCTAACACGTGCGCG", 
            "GTGGTGCATGGCCGTTCTT", "GTGGTGGTGCATGGCCGTTCTT", "TCGGTCGGGCTGGGGCGCGA")

new_names <- c("rsRNA-#122", "rsRNA-#146", "rsRNA-#216", "rsRNA-#217", "rsRNA-#250")  # 新名称

# 创建一个命名映射
name_mapping <- setNames(new_names, rsRNAs_to_annotate)

# 提取需要的 rsRNA 数据并设置顺序
heatmap_data <- merged_data %>%
  filter(rsRNA %in% rsRNAs_to_annotate) %>%
  select(rsRNA, log2control_median, log2treatment_median) %>%
  pivot_longer(cols = starts_with("log2"), names_to = "Condition", values_to = "Expression") %>%
  mutate(rsRNA = factor(name_mapping[rsRNA], levels = new_names))  # 使用新名称

# 定义颜色从白色到深蓝色
start_color <- rgb(255/255, 255/255, 255/255)  # 白色
end_color <- rgb(0/255, 51/255, 102/255)        # 深蓝色

# 生成渐变
num_colors <- 10
gradient_colors <- colorRampPalette(c(start_color, end_color))(num_colors)

# 绘制热图
heatmap_plot <- ggplot(heatmap_data, aes(x = rsRNA, y = Condition, fill = Expression)) +
  geom_tile(color = "white") +
  scale_fill_gradientn(colors = gradient_colors, 
                       name = "Expression Level") +
  labs(x = "rsRNA", y = "Condition", title = "Heatmap of rsRNA Expression") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 保存热图
ggsave("./MDD_hc_rsRNA/5rsRNA_heatmap.pdf", plot = heatmap_plot, device = "pdf", dpi = 600, width = 8, height = 6)
```
















## 1. other_mirna

2. 标注control_vs_treatment，手动标注
```r
setwd("D:/adult_dep/DE2")
df <- read.csv("./input/rpm_50%_log_mirna.csv", header = TRUE, row.names = 1)
df <- df[, df[1, ] %in% c(1, 2, 3)] # 一开始少了bp人群
colnames(df) <- ifelse(df[1, ] %in% c(2, 3), 
                       paste0("control_", colnames(df)), 
                       paste0("treatment_", colnames(df)))
control_columns <- colnames(df)[grepl("^control_", colnames(df))] 
treatment_columns <- colnames(df)[grepl("^treatment_", colnames(df))]        


# 删除第一行
df <- df[-1, ]
# 查看修改后的数据框
head(df)
write.csv(df, file = "./input_rename/other_mirna.csv")
```

3. 找到每个mirna的treatment_vs_control中位数，计算FC
```bash
cd /mnt/d/adult_dep/DE2
# 记得加“mirna”
python median.py ./input_rename/other_mirna.csv ./other_mirna/median.csv
```


② 缺失值替换为NA
```r
setwd("D:/adult_dep/DE2")
data <- read.csv("./input_rename/other_mirna.csv") 
data[data == 0] <- NA
write.csv(data, "./other_mirna/mirna_NA.csv", row.names = FALSE)  
```



③ Wilcoxon秩和检验（相关样本）  

秩和检验（Wilcoxon秩和检验）是一种用于比较两组相关或无关样本的非参数检验方法。它不要求数据符合正态分布，因此适用于那些无法满足正态性假设的情况。Wilcoxon秩和检验（相关样本）： 适用于配对设计，比如同一组个体在两个不同时间点的观测。在这种情况下，对每对相关样本进行秩和的比较。    
* 原理：在Wilcoxon秩和检验中，基本思想是将每个组的观测值按大小排序，然后为每个观测值分配一个秩次，最后通过比较秩和来判断两组样本是否来自同一分布。这个检验的零假设是两组样本来自相同的分布。  


```r
# 设置工作目录
setwd("D:/adult_dep/DE2")

# 读取数据
data <- read.csv("./other_mirna/mirna_NA.csv")

# 加载dplyr包
library(dplyr)

# 获取控制组和治疗组的列
control_columns <- grep("^control_", names(data), value = TRUE)
treatment_columns <- grep("^treatment_", names(data), value = TRUE)


# 执行 Wilcoxon 秩和检验
p_values_matrix <- matrix(NA, nrow = length(data$miRNA), ncol = 3)
colnames(p_values_matrix) <- c("miRNA", "p_value", "Q_value")

for (i in seq_along(data$miRNA)) {
  miRNA <- data$miRNA[i]
  control_values <- unlist(data[data$miRNA == miRNA, control_columns])
  treatment_values <- unlist(data[data$miRNA == miRNA, treatment_columns])
  
  # 执行 Wilcoxon 秩和检验
  result <- try(wilcox.test(treatment_values, control_values))

  # 存储 p-value
  if (!inherits(result, "try-error")) {
    p_values_matrix[i, 1] <- miRNA
    p_values_matrix[i, 2] <- result$p.value
  } else {
    cat("Error for miRNA:", miRNA, "\n")
  }
}

# 使用 Benjamini-Hochberg 方法进行多重比较校正
p_values_matrix[, 3] <- p.adjust(p_values_matrix[, 2], method = "BH")

# 输出结果
print(p_values_matrix)
# 将 p_values_matrix 保存到 CSV 文件
write.csv(p_values_matrix, file = "./other_mirna/de_values_matrix.csv", row.names = FALSE)
```

④ 合并文件
```r
library(dplyr)

median_data <- read.csv("./other_mirna/median.csv")
de_values_data <- read.csv("./other_mirna/de_values_matrix.csv")
# 基于 miRNA 列合并数据
merged_data <- inner_join(median_data, de_values_data, by = "miRNA")

# 查看合并后的数据
print(merged_data)
write.csv(merged_data, file = "./other_mirna/result.csv", row.names = FALSE)

# padj 小于 0.05 并且 Log2FC 大于 1（2倍） 或者小于 -1（1/2倍）
diff_gene <- subset(merged_data, Q_value < 0.05 & abs(log2fold_change) > 1)
write.csv(diff_gene, file="./other_mirna/sig_q0.05_fc2.csv", quote = F)
```

```r
library(dplyr)

median_data <- read.csv("./other_mirna/median.csv")
de_values_data <- read.csv("./other_mirna/de_values_matrix.csv")
# 基于 miRNA 列合并数据
merged_data <- inner_join(median_data, de_values_data, by = "miRNA")

# 查看合并后的数据
print(merged_data)
write.csv(merged_data, file = "./other_mirna/result.csv", row.names = FALSE)

diff_gene <- subset(merged_data, Q_value < 0.05 & abs(log2fold_change) > 0.58)
write.csv(diff_gene, file="./other_mirna/sig_q0.05_fc15.csv", quote = F)
```


⑤ 绘图
```r
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}
if (!requireNamespace("pheatmap", quietly = TRUE)) {
  install.packages("pheatmap")
}

# 加载包
library(ggplot2)
library(pheatmap)

merged_data$lgpadj <- -log10(merged_data$Q_value)

# 绘制火山图
volcano_plot <- ggplot(merged_data, aes(x = log2fold_change, y = lgpadj)) +
  geom_point(aes(color = ifelse(Q_value < 0.05 & log2fold_change > 1, "Up", 
                                ifelse(Q_value < 0.05 & log2fold_change < -1, "Down", "Not Sig"))), size = 2) +
  scale_color_manual(values = c("Up" = "#D32F2F", "Down" = "#1976D2", "Not Sig" = "grey70")) +
  labs(x = "Log2 Fold Change", y = "-Log10(padj)", title = "Volcano Plot") +
  theme_minimal(base_size = 14) +
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(),  # 去掉主网格线
        panel.grid.minor = element_blank(),  # 去掉次网格线
        axis.line = element_line(color = "black")) + # 显示横纵坐标轴
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", size = 0.5) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black", size = 0.5)

ggsave("./other_mirna/volcano_plot.pdf", plot = volcano_plot, device = "pdf", dpi = 600, width = 8, height = 6)
```

```r
# 加载包
library(ggplot2)
library(pheatmap)
setwd("D:/adult_dep/DE2")

merged_data <- read.csv("./other_mirna/result.csv")
merged_data$lgpadj <- -log10(merged_data$Q_value)

# 绘制火山图
volcano_plot <- ggplot(merged_data, aes(x = log2fold_change, y = lgpadj)) +
  geom_point(aes(color = ifelse(Q_value < 0.05 & log2fold_change > 0.58, "Up", 
                                ifelse(Q_value < 0.05 & log2fold_change < -0.58, "Down", "Not Sig"))), size = 2) +
  scale_color_manual(values = c("Up" = "#D32F2F", "Down" = "#1976D2", "Not Sig" = "grey70")) +
  labs(x = "Log2 Fold Change", y = "-Log10(padj)", title = "Volcano Plot") +
  theme_minimal(base_size = 14) +
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(),  # 去掉主网格线
        panel.grid.minor = element_blank(),  # 去掉次网格线
        axis.line = element_line(color = "black"), # 显示横纵坐标轴
        axis.ticks = element_line(color = "black"),  # 设置坐标轴刻度线颜色
        axis.ticks.length = unit(0.25, "cm")) +  # 设置坐标轴刻度线长度
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", size = 0.5) +
  geom_vline(xintercept = c(-0.58, 0.58), linetype = "dashed", color = "black", size = 0.5)

ggsave("./other_mirna/other_mirna_volcano_plot.pdf", plot = volcano_plot, device = "pdf", dpi = 600, width = 8, height = 6)
```

```r
# 图中标注名称
annotations <- data.frame(
  miRNA = c("hsa-mir-101-1_hsa-miR-101-3p", "hsa-mir-1-1_hsa-miR-1-3p", "hsa-mir-1-2_hsa-miR-1-3p","hsa-mir-127_hsa-miR-127-3p","hsa-mir-140_hsa-miR-140-3p","hsa-mir-15b_hsa-miR-15b-3p","hsa-mir-19b-1_hsa-miR-19b-3p","hsa-mir-20a_hsa-miR-20a-5p","hsa-mir-25_hsa-miR-25-3p","hsa-mir-4433b_hsa-miR-4433b-5p","hsa-mir-7-3_hsa-miR-7-5p"),
  new_label = c("miR-101-3p (mir-101-1)", "miR-1-3p (mir-1-1)", "miR-1-3p (mir-1-2)","miR-127-3p","miR-140-3p","miR-15b-3p","miR-19b-3p (mir-19b-1)","miR-20a-5p","miR-25-3p","miR-4433b-5p","miR-7-5p (mir-7-3)")
)

# 提取需要标注的 miRNA 数据
annotations_data <- merged_data[merged_data$miRNA %in% annotations$miRNA, ]

# 创建火山图
volcano_plot <- ggplot(merged_data, aes(x = log2fold_change, y = lgpadj)) +
  geom_point(aes(color = ifelse(Q_value < 0.05 & log2fold_change > 0.58, "Up", 
                                ifelse(Q_value < 0.05 & log2fold_change < -0.58, "Down", "Not Sig"))), size = 2) +
  scale_color_manual(values = c("Up" = "#D32F2F", "Down" = "#1976D2", "Not Sig" = "grey70")) +
  labs(x = "Log2 Fold Change", y = "-Log10(padj)", title = "Volcano Plot") +
  theme_minimal(base_size = 14) +
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", size = 0.5) +
  geom_vline(xintercept = c(-0.58, 0.58), linetype = "dashed", color = "black", size = 0.5) +
  geom_point(data = annotations_data, aes(color = "Highlighted"), size = 3, shape = 1, stroke = 1.5, fill = NA) +  # 添加圈圈
  geom_text(data = merge(annotations_data, annotations, by = "miRNA"), aes(label = new_label), 
            vjust = -1, hjust = 0.5, size = 3, color = "black")  # 使用新的标签


ggsave("./other_mirna/label_other_mirna_volcano_plot.pdf", plot = volcano_plot, device = "pdf", dpi = 600, width = 8, height = 6)

```

```r
# 绘制热图
library(ggplot2)
library(tidyr)
library(dplyr)

# 读取数据
merged_data <- read.csv("./other_mirna/result.csv")

# 指定要绘制热图的 miRNA
miRNAs_to_annotate <- c("hsa-mir-101-1_hsa-miR-101-3p", "hsa-mir-140_hsa-miR-140-3p", "hsa-mir-25_hsa-miR-25-3p",
                        "hsa-mir-20a_hsa-miR-20a-5p","hsa-mir-7-3_hsa-miR-7-5p","hsa-mir-15b_hsa-miR-15b-3p", "hsa-mir-19b-1_hsa-miR-19b-3p")

# 提取需要的 miRNA 数据并设置顺序
heatmap_data <- merged_data %>%
  filter(miRNA %in% miRNAs_to_annotate) %>%
  select(miRNA, log2control_median, log2treatment_median) %>%
  pivot_longer(cols = starts_with("log2"), names_to = "Condition", values_to = "Expression") %>%
  mutate(miRNA = factor(miRNA, levels = miRNAs_to_annotate))  # 设置因子顺序

# 定义颜色从白色到深蓝色
start_color <- rgb(255/255, 255/255, 255/255)  # 白色
end_color <- rgb(0/255, 51/255, 102/255)        # 深蓝色

# 生成渐变
num_colors <- 10
gradient_colors <- colorRampPalette(c(start_color, end_color))(num_colors)

# 绘制热图
heatmap_plot <- ggplot(heatmap_data, aes(x = miRNA, y = Condition, fill = Expression)) +
  geom_tile(color = "white") +
  scale_fill_gradientn(colors = gradient_colors, 
                       name = "Expression Level") +
  labs(x = "miRNA", y = "Condition", title = "Heatmap of miRNA Expression") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 保存热图
ggsave("./other_mirna/5mirna_heatmap.pdf", plot = heatmap_plot, device = "pdf", dpi = 600, width = 8, height = 6)
```



## 2. other_pirna

2. 标注control_vs_treatment，手动标注
```r
setwd("D:/adult_dep/DE2")
df <- read.csv("./input/rpm_50%_log_pirna.csv", header = TRUE, row.names = 1)
df <- df[, df[1, ] %in% c(1, 2, 3)]
colnames(df) <- ifelse(df[1, ] %in% c(2, 3), 
                       paste0("control_", colnames(df)), 
                       paste0("treatment_", colnames(df)))
control_columns <- colnames(df)[grepl("^control_", colnames(df))]        

# 删除第一行
df <- df[-1, ]
# 查看修改后的数据框
head(df)
write.csv(df, file = "./input_rename/other_pirna.csv")
```

3. 找到每个mirna的treatment_vs_control中位数，计算FC
```bash
cd /mnt/d/adult_dep/DE2
# 记得加“pirna”
python median.py ./input_rename/other_pirna.csv ./other_pirna/median.csv
```


② 缺失值替换为NA
```r
setwd("D:/adult_dep/DE2")
data <- read.csv("./input_rename/other_pirna.csv") 
data[data == 0] <- NA
write.csv(data, "./other_mirna/pirna_NA.csv", row.names = FALSE)  
```



③ Wilcoxon秩和检验（相关样本）  

秩和检验（Wilcoxon秩和检验）是一种用于比较两组相关或无关样本的非参数检验方法。它不要求数据符合正态分布，因此适用于那些无法满足正态性假设的情况。Wilcoxon秩和检验（相关样本）： 适用于配对设计，比如同一组个体在两个不同时间点的观测。在这种情况下，对每对相关样本进行秩和的比较。    
* 原理：在Wilcoxon秩和检验中，基本思想是将每个组的观测值按大小排序，然后为每个观测值分配一个秩次，最后通过比较秩和来判断两组样本是否来自同一分布。这个检验的零假设是两组样本来自相同的分布。  


```r
# 设置工作目录
setwd("D:/adult_dep/DE2")

# 读取数据
data <- read.csv("./other_mirna/pirna_NA.csv")

# 加载dplyr包
library(dplyr)

# 获取控制组和治疗组的列
control_columns <- grep("^control_", names(data), value = TRUE)
treatment_columns <- grep("^treatment_", names(data), value = TRUE)
# 835

# 执行 Wilcoxon 秩和检验
p_values_matrix <- matrix(NA, nrow = length(data$piRNA), ncol = 3)
colnames(p_values_matrix) <- c("piRNA", "p_value", "Q_value")

for (i in seq_along(data$piRNA)) {
  piRNA <- data$piRNA[i]
  control_values <- unlist(data[data$piRNA == piRNA, control_columns])
  treatment_values <- unlist(data[data$piRNA == piRNA, treatment_columns])
  
  # 执行 Wilcoxon 秩和检验
  result <- try(wilcox.test(treatment_values, control_values))

  # 存储 p-value
  if (!inherits(result, "try-error")) {
    p_values_matrix[i, 1] <- piRNA
    p_values_matrix[i, 2] <- result$p.value
  } else {
    cat("Error for piRNA:", piRNA, "\n")
  }
}

# 使用 Benjamini-Hochberg 方法进行多重比较校正
p_values_matrix[, 3] <- p.adjust(p_values_matrix[, 2], method = "BH")

# 输出结果
print(p_values_matrix)
# 将 p_values_matrix 保存到 CSV 文件
write.csv(p_values_matrix, file = "./other_pirna/de_values_matrix.csv", row.names = FALSE)
```

④ 合并文件
```r
library(dplyr)

median_data <- read.csv("./other_pirna/median.csv")
de_values_data <- read.csv("./other_pirna/de_values_matrix.csv")
# 基于 miRNA 列合并数据
merged_data <- inner_join(median_data, de_values_data, by = "piRNA")
# 查看合并后的数据
print(merged_data)
write.csv(merged_data, file = "./other_pirna/result.csv", row.names = FALSE)


# padj 小于 0.05 并且 Log2FC 大于 1（2倍） 或者小于 -1（1/2倍）
diff_gene <- subset(merged_data, Q_value < 0.05 & abs(log2fold_change) > 1)
write.csv(diff_gene, file="./other_pirna/sig_q0.05_fc2.csv", quote = F)
```

```r
library(dplyr)
setwd("D:/adult_dep/DE2")

merged_data <- read.csv("./other_pirna/result.csv")

diff_gene <- subset(merged_data, Q_value < 0.05 & abs(log2fold_change) > 0.58)
write.csv(diff_gene, file="./other_pirna/sig_q0.05_fc15.csv", quote = F)
```


⑤ 绘图
```r
# 加载包
library(ggplot2)
library(pheatmap)

merged_data$lgpadj <- -log10(merged_data$Q_value)

# 绘制火山图
volcano_plot <- ggplot(merged_data, aes(x = log2fold_change, y = lgpadj)) +
  geom_point(aes(color = ifelse(Q_value < 0.05 & log2fold_change > 1, "Up", 
                                ifelse(Q_value < 0.05 & log2fold_change < -1, "Down", "Not Sig"))), size = 2) +
  scale_color_manual(values = c("Up" = "#D32F2F", "Down" = "#1976D2", "Not Sig" = "grey70")) +
  labs(x = "Log2 Fold Change", y = "-Log10(padj)", title = "Volcano Plot") +
  theme_minimal(base_size = 14) +
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(),  # 去掉主网格线
        panel.grid.minor = element_blank(),  # 去掉次网格线
        axis.line = element_line(color = "black")) + # 显示横纵坐标轴
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", size = 0.5) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black", size = 0.5)

ggsave("./other_pirna/volcano_plot.pdf", plot = volcano_plot, device = "pdf", dpi = 600, width = 8, height = 6)
```




```r
# 加载包
library(ggplot2)
library(pheatmap)

merged_data <- read.csv("./other_pirna/result.csv")
merged_data$lgpadj <- -log10(merged_data$Q_value)

# 绘制火山图
volcano_plot <- ggplot(merged_data, aes(x = log2fold_change, y = lgpadj)) +
  geom_point(aes(color = ifelse(Q_value < 0.05 & log2fold_change > 0.58, "Up", 
                                ifelse(Q_value < 0.05 & log2fold_change < -0.58, "Down", "Not Sig"))), size = 2) +
  scale_color_manual(values = c("Up" = "#D32F2F", "Down" = "#1976D2", "Not Sig" = "grey70")) +
  labs(x = "Log2 Fold Change", y = "-Log10(padj)", title = "Volcano Plot") +
  theme_minimal(base_size = 14) +
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(),  # 去掉主网格线
        panel.grid.minor = element_blank(),  # 去掉次网格线
        axis.line = element_line(color = "black"), # 显示横纵坐标轴
        axis.ticks = element_line(color = "black"),  # 设置坐标轴刻度线颜色
        axis.ticks.length = unit(0.25, "cm")) +  # 设置坐标轴刻度线长度
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", size = 0.5) +
  geom_vline(xintercept = c(-0.58, 0.58), linetype = "dashed", color = "black", size = 0.5)

ggsave("./other_pirna/other_pirna_volcano_plot.pdf", plot = volcano_plot, device = "pdf", dpi = 600, width = 8, height = 6)
```
* 图中标注

```r
# 图中标注名称
miRNAs_to_annotate <- c("piR-hsa-168711", "piR-hsa-2482471","piR-hsa-334268","piR-hsa-35389","piR-hsa-3643462")

# 提取需要标注的 miRNA 数据
annotations <- merged_data[merged_data$piRNA %in% miRNAs_to_annotate, ]

# 创建火山图
volcano_plot <- ggplot(merged_data, aes(x = log2fold_change, y = lgpadj)) +
  geom_point(aes(color = ifelse(Q_value < 0.05 & log2fold_change > 0.58, "Up", 
                                ifelse(Q_value < 0.05 & log2fold_change < -0.58, "Down", "Not Sig"))), size = 2) +
  scale_color_manual(values = c("Up" = "#D32F2F", "Down" = "#1976D2", "Not Sig" = "grey70")) +
  labs(x = "Log2 Fold Change", y = "-Log10(padj)", title = "Volcano Plot") +
  theme_minimal(base_size = 14) +
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", size = 0.5) +
  geom_vline(xintercept = c(-0.58, 0.58), linetype = "dashed", color = "black", size = 0.5) +
  geom_point(data = annotations, aes(color = "Highlighted"), size = 3, shape = 1, stroke = 1.5, fill = NA) +  # 添加圈圈
  geom_text(data = annotations, aes(label = piRNA), vjust = -1, hjust = 0.5, size = 3, color = "black")


ggsave("./other_pirna/label_other_pirna_volcano_plot.pdf", plot = volcano_plot, device = "pdf", dpi = 600, width = 8, height = 6)

```

```r
# 绘制热图
library(ggplot2)
library(tidyr)
library(dplyr)

# 读取数据
merged_data <- read.csv("./other_pirna/result.csv")

# 指定要绘制热图的 miRNA
piRNAs_to_annotate <- c("piR-hsa-168711", "piR-hsa-2482471","piR-hsa-334268","piR-hsa-35389")

# 提取需要的 piRNA 数据并设置顺序
heatmap_data <- merged_data %>%
  filter(piRNA %in% piRNAs_to_annotate) %>%
  select(piRNA, log2control_median, log2treatment_median) %>%
  pivot_longer(cols = starts_with("log2"), names_to = "Condition", values_to = "Expression") %>%
  mutate(piRNA = factor(piRNA, levels = piRNAs_to_annotate))  # 设置因子顺序

# 定义颜色从白色到深蓝色
start_color <- rgb(255/255, 255/255, 255/255)  # 白色
end_color <- rgb(0/255, 51/255, 102/255)        # 深蓝色

# 生成渐变
num_colors <- 10
gradient_colors <- colorRampPalette(c(start_color, end_color))(num_colors)

# 绘制热图
heatmap_plot <- ggplot(heatmap_data, aes(x = piRNA, y = Condition, fill = Expression)) +
  geom_tile(color = "white") +
  scale_fill_gradientn(colors = gradient_colors, 
                       name = "Expression Level") +
  labs(x = "piRNA", y = "Condition", title = "Heatmap of piRNA Expression") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 保存热图
ggsave("./other_pirna/5piRNA_heatmap.pdf", plot = heatmap_plot, device = "pdf", dpi = 600, width = 8, height = 6)
```


























## 3. other_rsrna

2. 标注control_vs_treatment，手动标注
```r
setwd("D:/adult_dep/DE2")
df <- read.csv("./input/rpm_50%_log_rsrna.csv", header = TRUE, row.names = 1)
df <- df[, df[1, ] %in% c(1, 2, 3)]
colnames(df) <- ifelse(df[1, ] %in% c(2, 3), 
                       paste0("control_", colnames(df)), 
                       paste0("treatment_", colnames(df)))
control_columns <- colnames(df)[grepl("^control_", colnames(df))]        

# 删除第一行
df <- df[-1, ]
# 查看修改后的数据框
head(df)
write.csv(df, file = "./input_rename/other_rsrna.csv")
```

3. 找到每个mirna的treatment_vs_control中位数，计算FC
```bash
cd /mnt/d/adult_dep/DE2
# 记得加“rsrna”
python median.py ./input_rename/other_rsrna.csv ./other_rsrna/median.csv
```


② 缺失值替换为NA
```r
setwd("D:/adult_dep/DE2")
data <- read.csv("./input_rename/other_rsrna.csv") 
data[data == 0] <- NA
write.csv(data, "./other_rsrna/rsrna_NA.csv", row.names = FALSE)  
```



③ Wilcoxon秩和检验（相关样本）  

秩和检验（Wilcoxon秩和检验）是一种用于比较两组相关或无关样本的非参数检验方法。它不要求数据符合正态分布，因此适用于那些无法满足正态性假设的情况。Wilcoxon秩和检验（相关样本）： 适用于配对设计，比如同一组个体在两个不同时间点的观测。在这种情况下，对每对相关样本进行秩和的比较。    
* 原理：在Wilcoxon秩和检验中，基本思想是将每个组的观测值按大小排序，然后为每个观测值分配一个秩次，最后通过比较秩和来判断两组样本是否来自同一分布。这个检验的零假设是两组样本来自相同的分布。  


```r
# 设置工作目录
setwd("D:/adult_dep/DE2")

# 读取数据
data <- read.csv("./other_rsrna/rsrna_NA.csv")

# 加载dplyr包
library(dplyr)

# 获取控制组和治疗组的列
control_columns <- grep("^control_", names(data), value = TRUE)
treatment_columns <- grep("^treatment_", names(data), value = TRUE)

# 执行 Wilcoxon 秩和检验
p_values_matrix <- matrix(NA, nrow = length(data$rsRNA), ncol = 3)
colnames(p_values_matrix) <- c("rsRNA", "p_value", "Q_value")

for (i in seq_along(data$rsRNA)) {
  rsRNA <- data$rsRNA[i]
  control_values <- unlist(data[data$rsRNA == rsRNA, control_columns])
  treatment_values <- unlist(data[data$rsRNA == rsRNA, treatment_columns])
  
  # 执行 Wilcoxon 秩和检验
  result <- try(wilcox.test(treatment_values, control_values))

  # 存储 p-value
  if (!inherits(result, "try-error")) {
    p_values_matrix[i, 1] <- rsRNA
    p_values_matrix[i, 2] <- result$p.value
  } else {
    cat("Error for rsRNA:", rsRNA, "\n")
  }
}

# 使用 Benjamini-Hochberg 方法进行多重比较校正
p_values_matrix[, 3] <- p.adjust(p_values_matrix[, 2], method = "BH")

# 输出结果
print(p_values_matrix)
write.csv(p_values_matrix, file = "./other_rsrna/de_values_matrix.csv", row.names = FALSE)
```

④ 合并文件
```r
library(dplyr)

median_data <- read.csv("./other_rsrna/median.csv")
de_values_data <- read.csv("./other_rsrna/de_values_matrix.csv")
# 基于 miRNA 列合并数据
merged_data <- inner_join(median_data, de_values_data, by = "rsRNA")

# 查看合并后的数据
print(merged_data)
write.csv(merged_data, file = "./other_rsrna/result.csv", row.names = FALSE)

# padj 小于 0.05 并且 Log2FC 大于 1（2倍） 或者小于 -1（1/2倍）
diff_gene <- subset(merged_data, Q_value < 0.05 & abs(log2fold_change) > 1)
write.csv(diff_gene, file="./other_rsrna/sig_q0.05_fc2.csv", quote = F)
```


```r
library(dplyr)
setwd("D:/adult_dep/DE2")

merged_data <- read.csv("./other_rsrna/result.csv")

diff_gene <- subset(merged_data, Q_value < 0.05 & abs(log2fold_change) > 0.58)
write.csv(diff_gene, file="./other_rsrna/sig_q0.05_fc15.csv", quote = F)
```


⑤ 绘图
```r
# 加载包
library(ggplot2)
library(pheatmap)

merged_data$lgpadj <- -log10(merged_data$Q_value)

# 绘制火山图
volcano_plot <- ggplot(merged_data, aes(x = log2fold_change, y = lgpadj)) +
  geom_point(aes(color = ifelse(Q_value < 0.05 & log2fold_change > 1, "Up", 
                                ifelse(Q_value < 0.05 & log2fold_change < -1, "Down", "Not Sig"))), size = 2) +
  scale_color_manual(values = c("Up" = "#D32F2F", "Down" = "#1976D2", "Not Sig" = "grey70")) +
  labs(x = "Log2 Fold Change", y = "-Log10(padj)", title = "Volcano Plot") +
  theme_minimal(base_size = 14) +
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(),  # 去掉主网格线
        panel.grid.minor = element_blank(),  # 去掉次网格线
        axis.line = element_line(color = "black")) + # 显示横纵坐标轴
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", size = 0.5) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black", size = 0.5)

ggsave("./other_rsrna/volcano_plot.pdf", plot = volcano_plot, device = "pdf", dpi = 600, width = 8, height = 6)
```
```r
# 加载包
library(ggplot2)
library(pheatmap)

merged_data <- read.csv("./other_rsrna/result.csv")
merged_data$lgpadj <- -log10(merged_data$Q_value)

# 绘制火山图
volcano_plot <- ggplot(merged_data, aes(x = log2fold_change, y = lgpadj)) +
  geom_point(aes(color = ifelse(Q_value < 0.05 & log2fold_change > 0.58, "Up", 
                                ifelse(Q_value < 0.05 & log2fold_change < -0.58, "Down", "Not Sig"))), size = 2) +
  scale_color_manual(values = c("Up" = "#D32F2F", "Down" = "#1976D2", "Not Sig" = "grey70")) +
  labs(x = "Log2 Fold Change", y = "-Log10(padj)", title = "Volcano Plot") +
  theme_minimal(base_size = 14) +
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(),  # 去掉主网格线
        panel.grid.minor = element_blank(),  # 去掉次网格线
        axis.line = element_line(color = "black"), # 显示横纵坐标轴
        axis.ticks = element_line(color = "black"),  # 设置坐标轴刻度线颜色
        axis.ticks.length = unit(0.25, "cm")) +  # 设置坐标轴刻度线长度
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", size = 0.5) +
  geom_vline(xintercept = c(-0.58, 0.58), linetype = "dashed", color = "black", size = 0.5)

ggsave("./other_rsrna/other_rsrna_volcano_plot.pdf", plot = volcano_plot, device = "pdf", dpi = 600, width = 8, height = 6)
```
* 图中标注

```r
# 图中标注名称
annotations <- data.frame(
  rsRNA = c("ACGGACCAAGGAGTCTAACACGTGCGCG", "GGAGTCTAACACGTGCGCG"),
  new_label = c("rsRNA-#29", "rsRNA-#182")
)

# 提取需要标注的 miRNA 数据
annotations_data <- merged_data[merged_data$rsRNA %in% annotations$rsRNA, ]

# 创建火山图
volcano_plot <- ggplot(merged_data, aes(x = log2fold_change, y = lgpadj)) +
  geom_point(aes(color = ifelse(Q_value < 0.05 & log2fold_change > 0.58, "Up", 
                                ifelse(Q_value < 0.05 & log2fold_change < -0.58, "Down", "Not Sig"))), size = 2) +
  scale_color_manual(values = c("Up" = "#D32F2F", "Down" = "#1976D2", "Not Sig" = "grey70")) +
  labs(x = "Log2 Fold Change", y = "-Log10(padj)", title = "Volcano Plot") +
  theme_minimal(base_size = 14) +
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", size = 0.5) +
  geom_vline(xintercept = c(-0.58, 0.58), linetype = "dashed", color = "black", size = 0.5) +
  geom_point(data = annotations_data, aes(color = "Highlighted"), size = 3, shape = 1, stroke = 1.5, fill = NA) +  # 添加圈圈
  geom_text(data = merge(annotations_data, annotations, by = "rsRNA"), aes(label = new_label), 
            vjust = -1, hjust = 0.5, size = 3, color = "black")  # 使用新的标签


ggsave("./other_rsrna/label_other_rsrna_volcano_plot.pdf", plot = volcano_plot, device = "pdf", dpi = 600, width = 8, height = 6)

```














## 4. other_tsrna

2. 标注control_vs_treatment，手动标注
```r
setwd("D:/adult_dep/DE2")
df <- read.csv("./input/rpm_50%_log_tsrna.csv", header = TRUE, row.names = 1)
df <- df[, df[1, ] %in% c(1, 2, 3)]
colnames(df) <- ifelse(df[1, ] %in% c(2, 3), 
                       paste0("control_", colnames(df)), 
                       paste0("treatment_", colnames(df)))
control_columns <- colnames(df)[grepl("^control_", colnames(df))] 
treatment_columns <- colnames(df)[grepl("^treatment_", colnames(df))]        


# 删除第一行
df <- df[-1, ]
# 查看修改后的数据框
head(df)
write.csv(df, file = "./input_rename/other_tsrna.csv")
```

3. 找到每个mirna的treatment_vs_control中位数，计算FC
```bash
cd /mnt/d/adult_dep/DE2
# 记得加“tsrna”
python median.py ./input_rename/other_tsrna.csv ./other_tsrna/median.csv
```


② 缺失值替换为NA
```r
setwd("D:/adult_dep/DE2")
data <- read.csv("./input_rename/other_tsrna.csv") 
data[data == 0] <- NA
write.csv(data, "./other_tsrna/tsrna_NA.csv", row.names = FALSE)  
```



③ Wilcoxon秩和检验（相关样本）  

秩和检验（Wilcoxon秩和检验）是一种用于比较两组相关或无关样本的非参数检验方法。它不要求数据符合正态分布，因此适用于那些无法满足正态性假设的情况。Wilcoxon秩和检验（相关样本）： 适用于配对设计，比如同一组个体在两个不同时间点的观测。在这种情况下，对每对相关样本进行秩和的比较。    
* 原理：在Wilcoxon秩和检验中，基本思想是将每个组的观测值按大小排序，然后为每个观测值分配一个秩次，最后通过比较秩和来判断两组样本是否来自同一分布。这个检验的零假设是两组样本来自相同的分布。  


```r
# 设置工作目录
setwd("D:/adult_dep/DE2")

# 读取数据
data <- read.csv("./other_tsrna/tsrna_NA.csv")

# 加载dplyr包
library(dplyr)

# 获取控制组和治疗组的列
control_columns <- grep("^control_", names(data), value = TRUE)
treatment_columns <- grep("^treatment_", names(data), value = TRUE)

# 执行 Wilcoxon 秩和检验
p_values_matrix <- matrix(NA, nrow = length(data$tsRNA), ncol = 3)
colnames(p_values_matrix) <- c("tsRNA", "p_value", "Q_value")

for (i in seq_along(data$tsRNA)) {
  tsRNA <- data$tsRNA[i]
  control_values <- unlist(data[data$tsRNA == tsRNA, control_columns])
  treatment_values <- unlist(data[data$tsRNA == tsRNA, treatment_columns])
  
  # 执行 Wilcoxon 秩和检验
  result <- try(wilcox.test(treatment_values, control_values))

  # 存储 p-value
  if (!inherits(result, "try-error")) {
    p_values_matrix[i, 1] <- tsRNA
    p_values_matrix[i, 2] <- result$p.value
  } else {
    cat("Error for tsRNA:", tsRNA, "\n")
  }
}

# 使用 Benjamini-Hochberg 方法进行多重比较校正
p_values_matrix[, 3] <- p.adjust(p_values_matrix[, 2], method = "BH")

# 输出结果
print(p_values_matrix)
# 将 p_values_matrix 保存到 CSV 文件
write.csv(p_values_matrix, file = "./other_tsrna/de_values_matrix.csv", row.names = FALSE)
```

④ 合并文件
```r
library(dplyr)

median_data <- read.csv("./other_tsrna/median.csv")
de_values_data <- read.csv("./other_tsrna/de_values_matrix.csv")
# 基于 miRNA 列合并数据
merged_data <- inner_join(median_data, de_values_data, by = "tsRNA")

# 查看合并后的数据
print(merged_data)
write.csv(merged_data, file = "./other_tsrna/result.csv", row.names = FALSE)

# padj 小于 0.05 并且 Log2FC 大于 1（2倍） 或者小于 -1（1/2倍）
diff_gene <- subset(merged_data, Q_value < 0.05 & abs(log2fold_change) > 1)
write.csv(diff_gene, file="./other_tsrna/sig_q0.05_fc2.csv", quote = F)
```

```r
library(dplyr)
setwd("D:/adult_dep/DE2")

merged_data <- read.csv("./other_tsrna/result.csv")

diff_gene <- subset(merged_data, Q_value < 0.05 & abs(log2fold_change) > 0.58)
write.csv(diff_gene, file="./other_tsrna/sig_q0.05_fc15.csv", quote = F)
```


⑤ 绘图
```r
# 加载包
library(ggplot2)
library(pheatmap)

merged_data$lgpadj <- -log10(merged_data$Q_value)

# 绘制火山图
volcano_plot <- ggplot(merged_data, aes(x = log2fold_change, y = lgpadj)) +
  geom_point(aes(color = ifelse(Q_value < 0.05 & log2fold_change > 1, "Up", 
                                ifelse(Q_value < 0.05 & log2fold_change < -1, "Down", "Not Sig"))), size = 2) +
  scale_color_manual(values = c("Up" = "#D32F2F", "Down" = "#1976D2", "Not Sig" = "grey70")) +
  labs(x = "Log2 Fold Change", y = "-Log10(padj)", title = "Volcano Plot") +
  theme_minimal(base_size = 14) +
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(),  # 去掉主网格线
        panel.grid.minor = element_blank(),  # 去掉次网格线
        axis.line = element_line(color = "black")) + # 显示横纵坐标轴
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", size = 0.5) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black", size = 0.5)

ggsave("./other_tsrna/volcano_plot.pdf", plot = volcano_plot, device = "pdf", dpi = 600, width = 8, height = 6)
```
```r
# 加载包
library(ggplot2)
library(pheatmap)

merged_data <- read.csv("./other_tsrna/result.csv")
merged_data$lgpadj <- -log10(merged_data$Q_value)

# 绘制火山图
volcano_plot <- ggplot(merged_data, aes(x = log2fold_change, y = lgpadj)) +
  geom_point(aes(color = ifelse(Q_value < 0.05 & log2fold_change > 0.58, "Up", 
                                ifelse(Q_value < 0.05 & log2fold_change < -0.58, "Down", "Not Sig"))), size = 2) +
  scale_color_manual(values = c("Up" = "#D32F2F", "Down" = "#1976D2", "Not Sig" = "grey70")) +
  labs(x = "Log2 Fold Change", y = "-Log10(padj)", title = "Volcano Plot") +
  theme_minimal(base_size = 14) +
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(),  # 去掉主网格线
        panel.grid.minor = element_blank(),  # 去掉次网格线
        axis.line = element_line(color = "black"), # 显示横纵坐标轴
        axis.ticks = element_line(color = "black"),  # 设置坐标轴刻度线颜色
        axis.ticks.length = unit(0.25, "cm")) +  # 设置坐标轴刻度线长度
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", size = 0.5) +
  geom_vline(xintercept = c(-0.58, 0.58), linetype = "dashed", color = "black", size = 0.5)

ggsave("./other_tsrna/other_tsrna_volcano_plot.pdf", plot = volcano_plot, device = "pdf", dpi = 600, width = 8, height = 6)
```




* 图中标注

```r
# 图中标注名称
annotations <- data.frame(
  tsRNA = c("ATCCTGCCGACTACGCCA", "TCCCTGGTGGTCTAGTGGTTAGGATTCGGCGCT", "GCGGGAGACCGGGGTTCGATTCCCCGACGGGGAGC"),
  new_label = c("tsRNA-#17", "tsRNA-#83", "tsRNA-#51")
)

# 提取需要标注的 miRNA 数据
annotations_data <- merged_data[merged_data$tsRNA %in% annotations$tsRNA, ]

# 创建火山图
volcano_plot <- ggplot(merged_data, aes(x = log2fold_change, y = lgpadj)) +
  geom_point(aes(color = ifelse(Q_value < 0.05 & log2fold_change > 0.58, "Up", 
                                ifelse(Q_value < 0.05 & log2fold_change < -0.58, "Down", "Not Sig"))), size = 2) +
  scale_color_manual(values = c("Up" = "#D32F2F", "Down" = "#1976D2", "Not Sig" = "grey70")) +
  labs(x = "Log2 Fold Change", y = "-Log10(padj)", title = "Volcano Plot") +
  theme_minimal(base_size = 14) +
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", size = 0.5) +
  geom_vline(xintercept = c(-0.58, 0.58), linetype = "dashed", color = "black", size = 0.5) +
  geom_point(data = annotations_data, aes(color = "Highlighted"), size = 3, shape = 1, stroke = 1.5, fill = NA) +  # 添加圈圈
  geom_text(data = merge(annotations_data, annotations, by = "tsRNA"), aes(label = new_label), 
            vjust = -1, hjust = 0.5, size = 3, color = "black")  # 使用新的标签


ggsave("./other_tsrna/label_other_tsrna_volcano_plot.pdf", plot = volcano_plot, device = "pdf", dpi = 600, width = 8, height = 6)

```

```r
# 绘制热图
library(ggplot2)
library(tidyr)
library(dplyr)

# 读取数据
merged_data <- read.csv("./other_tsrna/result.csv")

# 指定需要标注的 tsRNA 和对应的新名称
tsRNAs_to_annotate <- c("ATCCTGCCGACTACGCCA", "TCCCTGGTGGTCTAGTGGTTAGGATTCGGCGCT", "GCGGGAGACCGGGGTTCGATTCCCCGACGGGGAGC")

new_names <- c("tsRNA-#17", "tsRNA-#83", "tsRNA-#51")  # 新名称

# 创建一个命名映射
name_mapping <- setNames(new_names, tsRNAs_to_annotate)

# 提取需要的 tsRNA 数据并设置顺序
heatmap_data <- merged_data %>%
  filter(tsRNA %in% tsRNAs_to_annotate) %>%
  select(tsRNA, log2control_median, log2treatment_median) %>%
  pivot_longer(cols = starts_with("log2"), names_to = "Condition", values_to = "Expression") %>%
  mutate(tsRNA = factor(name_mapping[tsRNA], levels = new_names))  # 使用新名称

# 定义颜色从白色到深蓝色
start_color <- rgb(255/255, 255/255, 255/255)  # 白色
end_color <- rgb(0/255, 51/255, 102/255)        # 深蓝色

# 生成渐变
num_colors <- 10
gradient_colors <- colorRampPalette(c(start_color, end_color))(num_colors)

# 绘制热图
heatmap_plot <- ggplot(heatmap_data, aes(x = tsRNA, y = Condition, fill = Expression)) +
  geom_tile(color = "white") +
  scale_fill_gradientn(colors = gradient_colors, 
                       name = "Expression Level") +
  labs(x = "tsRNA", y = "Condition", title = "Heatmap of tsRNA Expression") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 保存热图
ggsave("./other_tsrna/5tsRNA_heatmap.pdf", plot = heatmap_plot, device = "pdf", dpi = 600, width = 8, height = 6)
```






















### 1. bp_mirna

2. 标注control_vs_treatment，手动标注
```r
setwd("D:/adult_dep/DE2")
df <- read.csv("./input/rpm_50%_log_mirna.csv", header = TRUE, row.names = 1)
df <- df[, df[1, ] %in% c(1, 2)]
colnames(df) <- ifelse(df[1, ] == 2, 
                       paste0("control_", colnames(df)), 
                       paste0("treatment_", colnames(df)))
control_columns <- colnames(df)[grepl("^control_", colnames(df))]        

# 删除第一行
df <- df[-1, ]
# 查看修改后的数据框
head(df)
write.csv(df, file = "./input_rename/bp_mirna.csv")
```

3. 找到每个mirna的treatment_vs_control中位数，计算FC
```bash
cd /mnt/d/adult_dep/DE2
# 记得加“mirna”
python median.py ./input_rename/bp_mirna.csv ./bp_mirna/median.csv
python mean_mirna.py ./input_rename/bp_mirna.csv ./bp_mirna/mean.csv

```


② 缺失值替换为NA
```r
setwd("D:/adult_dep/DE2")
data <- read.csv("./input_rename/bp_mirna.csv") 
data[data == 0] <- NA
write.csv(data, "./bp_mirna/mirna_NA.csv", row.names = FALSE)  
```



③ Wilcoxon秩和检验（相关样本）  

秩和检验（Wilcoxon秩和检验）是一种用于比较两组相关或无关样本的非参数检验方法。它不要求数据符合正态分布，因此适用于那些无法满足正态性假设的情况。Wilcoxon秩和检验（相关样本）： 适用于配对设计，比如同一组个体在两个不同时间点的观测。在这种情况下，对每对相关样本进行秩和的比较。    
* 原理：在Wilcoxon秩和检验中，基本思想是将每个组的观测值按大小排序，然后为每个观测值分配一个秩次，最后通过比较秩和来判断两组样本是否来自同一分布。这个检验的零假设是两组样本来自相同的分布。  


```r
# 设置工作目录
setwd("D:/adult_dep/DE2")

# 读取数据
data <- read.csv("./bp_mirna/mirna_NA.csv")

# 加载dplyr包
library(dplyr)

# 获取控制组和治疗组的列
control_columns <- grep("^control_", names(data), value = TRUE)
treatment_columns <- grep("^treatment_", names(data), value = TRUE)


# 执行 Wilcoxon 秩和检验
p_values_matrix <- matrix(NA, nrow = length(data$miRNA), ncol = 3)
colnames(p_values_matrix) <- c("miRNA", "p_value", "Q_value")

for (i in seq_along(data$miRNA)) {
  miRNA <- data$miRNA[i]
  control_values <- unlist(data[data$miRNA == miRNA, control_columns])
  treatment_values <- unlist(data[data$miRNA == miRNA, treatment_columns])
  
  # 执行 Wilcoxon 秩和检验
  result <- try(wilcox.test(treatment_values, control_values))

  # 存储 p-value
  if (!inherits(result, "try-error")) {
    p_values_matrix[i, 1] <- miRNA
    p_values_matrix[i, 2] <- result$p.value
  } else {
    cat("Error for miRNA:", miRNA, "\n")
  }
}

# 使用 Benjamini-Hochberg 方法进行多重比较校正
p_values_matrix[, 3] <- p.adjust(p_values_matrix[, 2], method = "BH")

# 输出结果
print(p_values_matrix)
# 将 p_values_matrix 保存到 CSV 文件
write.csv(p_values_matrix, file = "./bp_mirna/de_values_matrix.csv", row.names = FALSE)
```

④ 合并文件
```r
library(dplyr)

median_data <- read.csv("./bp_mirna/median.csv")
de_values_data <- read.csv("./bp_mirna/de_values_matrix.csv")
# 基于 miRNA 列合并数据
merged_data <- inner_join(median_data, de_values_data, by = "miRNA")

# 查看合并后的数据
print(merged_data)
write.csv(merged_data, file = "./bp_mirna/result.csv", row.names = FALSE)

# padj 小于 0.05 并且 Log2FC 大于 1（2倍） 或者小于 -1（1/2倍）
diff_gene <- subset(merged_data, Q_value < 0.05 & abs(log2fold_change) > 1)
write.csv(diff_gene, file="./bp_mirna/sig_q0.05_fc2.csv", quote = F)
```

```r
library(dplyr)

median_data <- read.csv("./bp_mirna/median.csv")
de_values_data <- read.csv("./bp_mirna/de_values_matrix.csv")
# 基于 miRNA 列合并数据
merged_data <- inner_join(median_data, de_values_data, by = "miRNA")


diff_gene <- subset(merged_data, Q_value < 0.05 & abs(log2fold_change) > 0.58)
write.csv(diff_gene, file="./bp_mirna/sig_q0.05_fc15.csv", quote = F)



## mean
mean_data <- read.csv("./bp_mirna/mean.csv")
de_values_data <- read.csv("./bp_mirna/de_values_matrix.csv")
# 基于 miRNA 列合并数据
merged_data <- inner_join(mean_data, de_values_data, by = "miRNA")
diff_gene <- subset(merged_data, Q_value < 0.05 & abs(log2fold_change) > 0.58)
write.csv(diff_gene, file="./bp_mirna/mean_sig_q0.05_fc15.csv", quote = F)
```









⑤ 绘图
```r
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}
if (!requireNamespace("pheatmap", quietly = TRUE)) {
  install.packages("pheatmap")
}

# 加载包
library(ggplot2)
library(pheatmap)

merged_data$lgpadj <- -log10(merged_data$Q_value)

# 绘制火山图
volcano_plot <- ggplot(merged_data, aes(x = log2fold_change, y = lgpadj)) +
  geom_point(aes(color = ifelse(Q_value < 0.05 & log2fold_change > 1, "Up", 
                                ifelse(Q_value < 0.05 & log2fold_change < -1, "Down", "Not Sig"))), size = 2) +
  scale_color_manual(values = c("Up" = "#D32F2F", "Down" = "#1976D2", "Not Sig" = "grey70")) +
  labs(x = "Log2 Fold Change", y = "-Log10(padj)", title = "Volcano Plot") +
  theme_minimal(base_size = 14) +
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(),  # 去掉主网格线
        panel.grid.minor = element_blank(),  # 去掉次网格线
        axis.line = element_line(color = "black")) + # 显示横纵坐标轴
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", size = 0.5) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black", size = 0.5)

ggsave("./bp_mirna/volcano_plot.pdf", plot = volcano_plot, device = "pdf", dpi = 600, width = 8, height = 6)
```

```r
# 加载包
library(ggplot2)
library(pheatmap)
setwd("D:/adult_dep/DE2")

merged_data <- read.csv("./bp_mirna/result.csv")
merged_data$lgpadj <- -log10(merged_data$Q_value)

# 绘制火山图
volcano_plot <- ggplot(merged_data, aes(x = log2fold_change, y = lgpadj)) +
  geom_point(aes(color = ifelse(Q_value < 0.05 & log2fold_change > 0.58, "Up", 
                                ifelse(Q_value < 0.05 & log2fold_change < -0.58, "Down", "Not Sig"))), size = 2) +
  scale_color_manual(values = c("Up" = "#D32F2F", "Down" = "#1976D2", "Not Sig" = "grey70")) +
  labs(x = "Log2 Fold Change", y = "-Log10(padj)", title = "Volcano Plot") +
  theme_minimal(base_size = 14) +
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(),  # 去掉主网格线
        panel.grid.minor = element_blank(),  # 去掉次网格线
        axis.line = element_line(color = "black"), # 显示横纵坐标轴
        axis.ticks = element_line(color = "black"),  # 设置坐标轴刻度线颜色
        axis.ticks.length = unit(0.25, "cm")) +  # 设置坐标轴刻度线长度
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", size = 0.5) +
  geom_vline(xintercept = c(-0.58, 0.58), linetype = "dashed", color = "black", size = 0.5)

ggsave("./bp_mirna/bp_mirna_volcano_plot.pdf", plot = volcano_plot, device = "pdf", dpi = 600, width = 8, height = 6)
```

```r
# 图中标注名称
miRNAs_to_annotate <- c("hsa-mir-20a_hsa-miR-20a-5p", "hsa-mir-410_hsa-miR-410-3p", 
                         "hsa-mir-4433b_hsa-miR-4433b-5p", "hsa-mir-93_hsa-miR-93-5p", "hsa-mir-92a-1_hsa-miR-92a-1-5p","hsa-mir-532_hsa-miR-532-5p","hsa-mir-1827_hsa-miR-1827")

# 提取需要标注的 miRNA 数据
annotations <- merged_data[merged_data$miRNA %in% miRNAs_to_annotate, ]

# 创建火山图
volcano_plot <- ggplot(merged_data, aes(x = log2fold_change, y = lgpadj)) +
  geom_point(aes(color = ifelse(Q_value < 0.05 & log2fold_change > 0.58, "Up", 
                                ifelse(Q_value < 0.05 & log2fold_change < -0.58, "Down", "Not Sig"))), size = 2) +
  scale_color_manual(values = c("Up" = "#D32F2F", "Down" = "#1976D2", "Not Sig" = "grey70")) +
  labs(x = "Log2 Fold Change", y = "-Log10(padj)", title = "Volcano Plot") +
  theme_minimal(base_size = 14) +
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", size = 0.5) +
  geom_vline(xintercept = c(-0.58, 0.58), linetype = "dashed", color = "black", size = 0.5) +
  geom_point(data = annotations, aes(color = "Highlighted"), size = 3, shape = 1, stroke = 1.5, fill = NA) +  # 添加圈圈
  geom_text(data = annotations, aes(label = miRNA), vjust = -1, hjust = 0.5, size = 3, color = "black")


ggsave("./bp_mirna/label_bp_mirna_volcano_plot.pdf", plot = volcano_plot, device = "pdf", dpi = 600, width = 8, height = 6)
ggsave("./bp_mirna/label_bp_mirna_volcano_plot2.pdf", plot = volcano_plot, device = "pdf", dpi = 600, width = 16, height = 12)


```

```r
# 绘制热图
library(ggplot2)
library(tidyr)
library(dplyr)

# 读取数据
merged_data <- read.csv("./bp_mirna/result.csv")

# 指定要绘制热图的 miRNA
miRNAs_to_annotate <- c("hsa-mir-28_hsa-miR-28-3p", "hsa-mir-410_hsa-miR-410-3p", 
                         "hsa-mir-4433b_hsa-miR-4433b-5p")

# 提取需要的 miRNA 数据并设置顺序
heatmap_data <- merged_data %>%
  filter(miRNA %in% miRNAs_to_annotate) %>%
  select(miRNA, log2control_median, log2treatment_median) %>%
  pivot_longer(cols = starts_with("log2"), names_to = "Condition", values_to = "Expression") %>%
  mutate(miRNA = factor(miRNA, levels = miRNAs_to_annotate))  # 设置因子顺序

# 定义颜色从白色到深蓝色
start_color <- rgb(255/255, 255/255, 255/255)  # 白色
end_color <- rgb(0/255, 51/255, 102/255)        # 深蓝色

# 生成渐变
num_colors <- 10
gradient_colors <- colorRampPalette(c(start_color, end_color))(num_colors)

# 绘制热图
heatmap_plot <- ggplot(heatmap_data, aes(x = miRNA, y = Condition, fill = Expression)) +
  geom_tile(color = "white") +
  scale_fill_gradientn(colors = gradient_colors, 
                       name = "Expression Level") +
  labs(x = "miRNA", y = "Condition", title = "Heatmap of miRNA Expression") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 保存热图
ggsave("./bp_mirna/5mirna_heatmap.pdf", plot = heatmap_plot, device = "pdf", dpi = 600, width = 8, height = 6)
```
























### 2. bp_pirna

2. 标注control_vs_treatment，手动标注
```r
setwd("D:/adult_dep/DE2")
df <- read.csv("./input/rpm_50%_log_pirna.csv", header = TRUE, row.names = 1)
df <- df[, df[1, ] %in% c(1, 2)]
colnames(df) <- ifelse(df[1, ] == 2, 
                       paste0("control_", colnames(df)), 
                       paste0("treatment_", colnames(df)))
control_columns <- colnames(df)[grepl("^control_", colnames(df))]        

# 删除第一行
df <- df[-1, ]
# 查看修改后的数据框
head(df)
write.csv(df, file = "./input_rename/bp_pirna.csv")
```

3. 找到每个mirna的treatment_vs_control中位数，计算FC
```bash
cd /mnt/d/adult_dep/DE2
# 记得加“pirna”
python median.py ./input_rename/bp_pirna.csv ./bp_pirna/median.csv
```


② 缺失值替换为NA
```r
setwd("D:/adult_dep/DE2")
data <- read.csv("./input_rename/bp_pirna.csv") 
data[data == 0] <- NA
write.csv(data, "./bp_pirna/pirna_NA.csv", row.names = FALSE)  
```



③ Wilcoxon秩和检验（相关样本）  

秩和检验（Wilcoxon秩和检验）是一种用于比较两组相关或无关样本的非参数检验方法。它不要求数据符合正态分布，因此适用于那些无法满足正态性假设的情况。Wilcoxon秩和检验（相关样本）： 适用于配对设计，比如同一组个体在两个不同时间点的观测。在这种情况下，对每对相关样本进行秩和的比较。    
* 原理：在Wilcoxon秩和检验中，基本思想是将每个组的观测值按大小排序，然后为每个观测值分配一个秩次，最后通过比较秩和来判断两组样本是否来自同一分布。这个检验的零假设是两组样本来自相同的分布。  


```r
# 设置工作目录
setwd("D:/adult_dep/DE2")

# 读取数据
data <- read.csv("./bp_pirna/pirna_NA.csv")

# 加载dplyr包
library(dplyr)

# 获取控制组和治疗组的列
control_columns <- grep("^control_", names(data), value = TRUE)
treatment_columns <- grep("^treatment_", names(data), value = TRUE)
# 835

# 执行 Wilcoxon 秩和检验
p_values_matrix <- matrix(NA, nrow = length(data$piRNA), ncol = 3)
colnames(p_values_matrix) <- c("piRNA", "p_value", "Q_value")

for (i in seq_along(data$piRNA)) {
  piRNA <- data$piRNA[i]
  control_values <- unlist(data[data$piRNA == piRNA, control_columns])
  treatment_values <- unlist(data[data$piRNA == piRNA, treatment_columns])
  
  # 执行 Wilcoxon 秩和检验
  result <- try(wilcox.test(treatment_values, control_values))

  # 存储 p-value
  if (!inherits(result, "try-error")) {
    p_values_matrix[i, 1] <- piRNA
    p_values_matrix[i, 2] <- result$p.value
  } else {
    cat("Error for piRNA:", piRNA, "\n")
  }
}

# 使用 Benjamini-Hochberg 方法进行多重比较校正
p_values_matrix[, 3] <- p.adjust(p_values_matrix[, 2], method = "BH")

# 输出结果
print(p_values_matrix)
# 将 p_values_matrix 保存到 CSV 文件
write.csv(p_values_matrix, file = "./bp_pirna/de_values_matrix.csv", row.names = FALSE)
```

④ 合并文件
```r
library(dplyr)

median_data <- read.csv("./bp_pirna/median.csv")
de_values_data <- read.csv("./bp_pirna/de_values_matrix.csv")
# 基于 miRNA 列合并数据
merged_data <- inner_join(median_data, de_values_data, by = "piRNA")

# 查看合并后的数据
print(merged_data)
write.csv(merged_data, file = "./bp_pirna/result.csv", row.names = FALSE)

# padj 小于 0.05 并且 Log2FC 大于 1（2倍） 或者小于 -1（1/2倍）
diff_gene <- subset(merged_data, Q_value < 0.05 & abs(log2fold_change) > 1)
write.csv(diff_gene, file="./bp_pirna/sig_q0.05_fc2.csv", quote = F)
```

```r
library(dplyr)
setwd("D:/adult_dep/DE2")

merged_data <- read.csv("./bp_pirna/result.csv")

diff_gene <- subset(merged_data, Q_value < 0.05 & abs(log2fold_change) > 0.58)
write.csv(diff_gene, file="./bp_pirna/sig_q0.05_fc15.csv", quote = F)
```




⑤ 绘图
```r
# 加载包
library(ggplot2)
library(pheatmap)

merged_data$lgpadj <- -log10(merged_data$Q_value)

# 绘制火山图
volcano_plot <- ggplot(merged_data, aes(x = log2fold_change, y = lgpadj)) +
  geom_point(aes(color = ifelse(Q_value < 0.05 & log2fold_change > 1, "Up", 
                                ifelse(Q_value < 0.05 & log2fold_change < -1, "Down", "Not Sig"))), size = 2) +
  scale_color_manual(values = c("Up" = "#D32F2F", "Down" = "#1976D2", "Not Sig" = "grey70")) +
  labs(x = "Log2 Fold Change", y = "-Log10(padj)", title = "Volcano Plot") +
  theme_minimal(base_size = 14) +
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(),  # 去掉主网格线
        panel.grid.minor = element_blank(),  # 去掉次网格线
        axis.line = element_line(color = "black")) + # 显示横纵坐标轴
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", size = 0.5) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black", size = 0.5)

ggsave("./bp_pirna/volcano_plot.pdf", plot = volcano_plot, device = "pdf", dpi = 600, width = 8, height = 6)
```




```r
# 加载包
library(ggplot2)
library(pheatmap)

merged_data <- read.csv("./bp_pirna/result.csv")
merged_data$lgpadj <- -log10(merged_data$Q_value)

# 绘制火山图
volcano_plot <- ggplot(merged_data, aes(x = log2fold_change, y = lgpadj)) +
  geom_point(aes(color = ifelse(Q_value < 0.05 & log2fold_change > 0.58, "Up", 
                                ifelse(Q_value < 0.05 & log2fold_change < -0.58, "Down", "Not Sig"))), size = 2) +
  scale_color_manual(values = c("Up" = "#D32F2F", "Down" = "#1976D2", "Not Sig" = "grey70")) +
  labs(x = "Log2 Fold Change", y = "-Log10(padj)", title = "Volcano Plot") +
  theme_minimal(base_size = 14) +
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(),  # 去掉主网格线
        panel.grid.minor = element_blank(),  # 去掉次网格线
        axis.line = element_line(color = "black"), # 显示横纵坐标轴
        axis.ticks = element_line(color = "black"),  # 设置坐标轴刻度线颜色
        axis.ticks.length = unit(0.25, "cm")) +  # 设置坐标轴刻度线长度
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", size = 0.5) +
  geom_vline(xintercept = c(-0.58, 0.58), linetype = "dashed", color = "black", size = 0.5)

ggsave("./bp_pirna/bp_pirna_volcano_plot.pdf", plot = volcano_plot, device = "pdf", dpi = 600, width = 8, height = 6)
```
* 图中标注

```r
# 图中标注名称
miRNAs_to_annotate <- c("piR-hsa-135526", "piR-hsa-147698","piR-hsa-151782","piR-hsa-184607","piR-hsa-789156")

# 提取需要标注的 miRNA 数据
annotations <- merged_data[merged_data$piRNA %in% miRNAs_to_annotate, ]

# 创建火山图
volcano_plot <- ggplot(merged_data, aes(x = log2fold_change, y = lgpadj)) +
  geom_point(aes(color = ifelse(Q_value < 0.05 & log2fold_change > 0.58, "Up", 
                                ifelse(Q_value < 0.05 & log2fold_change < -0.58, "Down", "Not Sig"))), size = 2) +
  scale_color_manual(values = c("Up" = "#D32F2F", "Down" = "#1976D2", "Not Sig" = "grey70")) +
  labs(x = "Log2 Fold Change", y = "-Log10(padj)", title = "Volcano Plot") +
  theme_minimal(base_size = 14) +
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", size = 0.5) +
  geom_vline(xintercept = c(-0.58, 0.58), linetype = "dashed", color = "black", size = 0.5) +
  geom_point(data = annotations, aes(color = "Highlighted"), size = 3, shape = 1, stroke = 1.5, fill = NA) +  # 添加圈圈
  geom_text(data = annotations, aes(label = piRNA), vjust = -1, hjust = 0.5, size = 3, color = "black")


ggsave("./bp_pirna/label_bp_pirna_volcano_plot.pdf", plot = volcano_plot, device = "pdf", dpi = 600, width = 8, height = 6)










###################### 第二次绘图  ######################
# 图中标注名称
miRNAs_to_annotate <- c("piR-hsa-35389")

# Extract the data for miRNAs to annotate
annotations <- merged_data[merged_data$piRNA %in% miRNAs_to_annotate, ]

# Add annotations for the top 5 miRNAs with the highest p-value (lowest significance)
top_p_values <- merged_data[order(merged_data$lgpadj, decreasing = TRUE), ]
top_p_annotations <- head(top_p_values, 3)

# Combine the annotations: predefined list and the top 5 p-value miRNAs
all_annotations <- rbind(annotations, top_p_annotations)

# Create volcano plot
volcano_plot <- ggplot(merged_data, aes(x = log2fold_change, y = lgpadj)) +
  geom_point(aes(color = ifelse(Q_value < 0.05 & log2fold_change > 0.58, "Up", 
                                ifelse(Q_value < 0.05 & log2fold_change < -0.58, "Down", "Not Sig"))), size = 2) +
  scale_color_manual(values = c("Up" = "#D32F2F", "Down" = "#1976D2", "Not Sig" = "grey70")) +
  labs(x = "Log2 Fold Change", y = "-Log10(padj)", title = "Volcano Plot") +
  theme_minimal(base_size = 14) +
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", size = 0.5) +
  geom_vline(xintercept = c(-0.58, 0.58), linetype = "dashed", color = "black", size = 0.5) +
  geom_point(data = annotations, aes(color = "Highlighted"), size = 3, shape = 1, stroke = 1.5, fill = NA) +  # Highlight miRNAs from predefined list
  geom_text(data = all_annotations, aes(label = piRNA), vjust = -1, hjust = 0.5, size = 3, color = "black")  # Add text labels for both predefined and top p-value miRNAs

# Display the volcano plot
print(volcano_plot)


ggsave("./bp_pirna/label_pirna_volcano_plot3.pdf", plot = volcano_plot, device = "pdf", dpi = 600, width = 16, height = 12)

ggsave("./bp_pirna/label_pirna_volcano_plot2.pdf", plot = volcano_plot, device = "pdf", dpi = 600, width = 8, height = 6)
```

```r
# 绘制热图
library(ggplot2)
library(tidyr)
library(dplyr)

# 读取数据
merged_data <- read.csv("./bp_pirna/result.csv")

# 指定要绘制热图的 miRNA
piRNAs_to_annotate <- c("piR-hsa-135526", "piR-hsa-147698","piR-hsa-151782","piR-hsa-184607","piR-hsa-789156")

# 提取需要的 piRNA 数据并设置顺序
heatmap_data <- merged_data %>%
  filter(piRNA %in% piRNAs_to_annotate) %>%
  select(piRNA, log2control_median, log2treatment_median) %>%
  pivot_longer(cols = starts_with("log2"), names_to = "Condition", values_to = "Expression") %>%
  mutate(piRNA = factor(piRNA, levels = piRNAs_to_annotate))  # 设置因子顺序

# 定义颜色从白色到深蓝色
start_color <- rgb(255/255, 255/255, 255/255)  # 白色
end_color <- rgb(0/255, 51/255, 102/255)        # 深蓝色

# 生成渐变
num_colors <- 10
gradient_colors <- colorRampPalette(c(start_color, end_color))(num_colors)

# 绘制热图
heatmap_plot <- ggplot(heatmap_data, aes(x = piRNA, y = Condition, fill = Expression)) +
  geom_tile(color = "white") +
  scale_fill_gradientn(colors = gradient_colors, 
                       name = "Expression Level") +
  labs(x = "piRNA", y = "Condition", title = "Heatmap of piRNA Expression") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 保存热图
ggsave("./bp_pirna/5piRNA_heatmap.pdf", plot = heatmap_plot, device = "pdf", dpi = 600, width = 8, height = 6)
```






## 3. bp_rsrna

2. 标注control_vs_treatment，手动标注
```r
setwd("D:/adult_dep/DE2")
df <- read.csv("./input/rpm_50%_log_rsrna.csv", header = TRUE, row.names = 1)
df <- df[, df[1, ] %in% c(1, 2)]
colnames(df) <- ifelse(df[1, ] == 2, 
                       paste0("control_", colnames(df)), 
                       paste0("treatment_", colnames(df)))
control_columns <- colnames(df)[grepl("^control_", colnames(df))]        

# 删除第一行
df <- df[-1, ]
# 查看修改后的数据框
head(df)
write.csv(df, file = "./input_rename/bp_rsrna.csv")
```

3. 找到每个mirna的treatment_vs_control中位数，计算FC
```bash
cd /mnt/d/adult_dep/DE2
# 记得加“rsrna”
python median.py ./input_rename/bp_rsrna.csv ./bp_rsrna/median.csv
```


② 缺失值替换为NA
```r
setwd("D:/adult_dep/DE2")
data <- read.csv("./input_rename/bp_rsrna.csv") 
data[data == 0] <- NA
write.csv(data, "./bp_rsrna/rsrna_NA.csv", row.names = FALSE)  
```



③ Wilcoxon秩和检验（相关样本）  

秩和检验（Wilcoxon秩和检验）是一种用于比较两组相关或无关样本的非参数检验方法。它不要求数据符合正态分布，因此适用于那些无法满足正态性假设的情况。Wilcoxon秩和检验（相关样本）： 适用于配对设计，比如同一组个体在两个不同时间点的观测。在这种情况下，对每对相关样本进行秩和的比较。    
* 原理：在Wilcoxon秩和检验中，基本思想是将每个组的观测值按大小排序，然后为每个观测值分配一个秩次，最后通过比较秩和来判断两组样本是否来自同一分布。这个检验的零假设是两组样本来自相同的分布。  


```r
# 设置工作目录
setwd("D:/adult_dep/DE2")

# 读取数据
data <- read.csv("./bp_rsrna/rsrna_NA.csv")

# 加载dplyr包
library(dplyr)

# 获取控制组和治疗组的列
control_columns <- grep("^control_", names(data), value = TRUE)
treatment_columns <- grep("^treatment_", names(data), value = TRUE)

# 执行 Wilcoxon 秩和检验
p_values_matrix <- matrix(NA, nrow = length(data$rsRNA), ncol = 3)
colnames(p_values_matrix) <- c("rsRNA", "p_value", "Q_value")

for (i in seq_along(data$rsRNA)) {
  rsRNA <- data$rsRNA[i]
  control_values <- unlist(data[data$rsRNA == rsRNA, control_columns])
  treatment_values <- unlist(data[data$rsRNA == rsRNA, treatment_columns])
  
  # 执行 Wilcoxon 秩和检验
  result <- try(wilcox.test(treatment_values, control_values))

  # 存储 p-value
  if (!inherits(result, "try-error")) {
    p_values_matrix[i, 1] <- rsRNA
    p_values_matrix[i, 2] <- result$p.value
  } else {
    cat("Error for rsRNA:", rsRNA, "\n")
  }
}

# 使用 Benjamini-Hochberg 方法进行多重比较校正
p_values_matrix[, 3] <- p.adjust(p_values_matrix[, 2], method = "BH")

# 输出结果
print(p_values_matrix)
write.csv(p_values_matrix, file = "./bp_rsrna/de_values_matrix.csv", row.names = FALSE)
```

④ 合并文件
```r
library(dplyr)

median_data <- read.csv("./bp_rsrna/median.csv")
de_values_data <- read.csv("./bp_rsrna/de_values_matrix.csv")
# 基于 miRNA 列合并数据
merged_data <- inner_join(median_data, de_values_data, by = "rsRNA")

# 查看合并后的数据
print(merged_data)
write.csv(merged_data, file = "./bp_rsrna/result.csv", row.names = FALSE)

# padj 小于 0.05 并且 Log2FC 大于 1（2倍） 或者小于 -1（1/2倍）
diff_gene <- subset(merged_data, Q_value < 0.05 & abs(log2fold_change) > 1)
write.csv(diff_gene, file="./bp_rsrna/sig_q0.05_fc2.csv", quote = F)
```



```r
library(dplyr)
setwd("D:/adult_dep/DE2")

merged_data <- read.csv("./bp_rsrna/result.csv")

diff_gene <- subset(merged_data, Q_value < 0.05 & abs(log2fold_change) > 0.58)
write.csv(diff_gene, file="./bp_rsrna/sig_q0.05_fc15.csv", quote = F)
```



⑤ 绘图
```r
# 加载包
library(ggplot2)
library(pheatmap)

merged_data$lgpadj <- -log10(merged_data$Q_value)

# 绘制火山图
volcano_plot <- ggplot(merged_data, aes(x = log2fold_change, y = lgpadj)) +
  geom_point(aes(color = ifelse(Q_value < 0.05 & log2fold_change > 1, "Up", 
                                ifelse(Q_value < 0.05 & log2fold_change < -1, "Down", "Not Sig"))), size = 2) +
  scale_color_manual(values = c("Up" = "#D32F2F", "Down" = "#1976D2", "Not Sig" = "grey70")) +
  labs(x = "Log2 Fold Change", y = "-Log10(padj)", title = "Volcano Plot") +
  theme_minimal(base_size = 14) +
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(),  # 去掉主网格线
        panel.grid.minor = element_blank(),  # 去掉次网格线
        axis.line = element_line(color = "black")) + # 显示横纵坐标轴
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", size = 0.5) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black", size = 0.5)

ggsave("./bp_rsrna/volcano_plot.pdf", plot = volcano_plot, device = "pdf", dpi = 600, width = 8, height = 6)
```

```r
# 加载包
library(ggplot2)
library(pheatmap)

merged_data <- read.csv("./bp_rsrna/result.csv")
merged_data$lgpadj <- -log10(merged_data$Q_value)

# 绘制火山图
volcano_plot <- ggplot(merged_data, aes(x = log2fold_change, y = lgpadj)) +
  geom_point(aes(color = ifelse(Q_value < 0.05 & log2fold_change > 0.58, "Up", 
                                ifelse(Q_value < 0.05 & log2fold_change < -0.58, "Down", "Not Sig"))), size = 2) +
  scale_color_manual(values = c("Up" = "#D32F2F", "Down" = "#1976D2", "Not Sig" = "grey70")) +
  labs(x = "Log2 Fold Change", y = "-Log10(padj)", title = "Volcano Plot") +
  theme_minimal(base_size = 14) +
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(),  # 去掉主网格线
        panel.grid.minor = element_blank(),  # 去掉次网格线
        axis.line = element_line(color = "black"), # 显示横纵坐标轴
        axis.ticks = element_line(color = "black"),  # 设置坐标轴刻度线颜色
        axis.ticks.length = unit(0.25, "cm")) +  # 设置坐标轴刻度线长度
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", size = 0.5) +
  geom_vline(xintercept = c(-0.58, 0.58), linetype = "dashed", color = "black", size = 0.5)

ggsave("./bp_rsrna/bp_rsrna_volcano_plot.pdf", plot = volcano_plot, device = "pdf", dpi = 600, width = 8, height = 6)
```




* 图中标注

```r
# 图中标注名称
annotations <- data.frame(
  rsRNA = c("ACAATGTAGGTAAGGGAAGTCG", "ACGGACCAAGGAGTCTAACACGTGCGCG", "CGACTCTTAGCGGTGGATCACTCGGCTCGTGCG","GGAACAATGTAGGTAAGGGA","GGAGTCTAACACGTGCGCG","TGGGAGACCGCCTGGGAATACCGGGTGCTGTAGGCTT", "GACTCTTAGCGGTGGATCACTCGGCTCGTG"),
  new_label = c("rsRNA-#18", "rsRNA-#29", "rsRNA-#95","rsRNA-#176","rsRNA-#182","rsRNA-#261", "rsRNA-#155")
)

# 提取需要标注的 miRNA 数据
annotations_data <- merged_data[merged_data$rsRNA %in% annotations$rsRNA, ]

# 创建火山图
volcano_plot <- ggplot(merged_data, aes(x = log2fold_change, y = lgpadj)) +
  geom_point(aes(color = ifelse(Q_value < 0.05 & log2fold_change > 0.58, "Up", 
                                ifelse(Q_value < 0.05 & log2fold_change < -0.58, "Down", "Not Sig"))), size = 2) +
  scale_color_manual(values = c("Up" = "#D32F2F", "Down" = "#1976D2", "Not Sig" = "grey70")) +
  labs(x = "Log2 Fold Change", y = "-Log10(padj)", title = "Volcano Plot") +
  theme_minimal(base_size = 14) +
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", size = 0.5) +
  geom_vline(xintercept = c(-0.58, 0.58), linetype = "dashed", color = "black", size = 0.5) +
  geom_point(data = annotations_data, aes(color = "Highlighted"), size = 3, shape = 1, stroke = 1.5, fill = NA) +  # 添加圈圈
  geom_text(data = merge(annotations_data, annotations, by = "rsRNA"), aes(label = new_label), 
            vjust = -1, hjust = 0.5, size = 3, color = "black")  # 使用新的标签


ggsave("./bp_rsrna/label_bp_rsrna_volcano_plot.pdf", plot = volcano_plot, device = "pdf", dpi = 600, width = 8, height = 6)

```

```r
# 绘制热图
library(ggplot2)
library(tidyr)
library(dplyr)

# 读取数据
merged_data <- read.csv("./bp_rsrna/result.csv")

# 指定需要标注的 tsRNA 和对应的新名称
rsRNAs_to_annotate <- c("ACAATGTAGGTAAGGGAAGTCG", "ACGGACCAAGGAGTCTAACACGTGCGCG", "CGACTCTTAGCGGTGGATCACTCGGCTCGTGCG","GGAACAATGTAGGTAAGGGA","GGAGTCTAACACGTGCGCG","TGGGAGACCGCCTGGGAATACCGGGTGCTGTAGGCTT", "GACTCTTAGCGGTGGATCACTCGGCTCGTG")

new_names <- c("rsRNA-#18", "rsRNA-#29", "rsRNA-#95","rsRNA-#176","rsRNA-#182","rsRNA-#261", "rsRNA-#155")

# 创建一个命名映射
name_mapping <- setNames(new_names, rsRNAs_to_annotate)

# 提取需要的 rsRNA 数据并设置顺序
heatmap_data <- merged_data %>%
  filter(rsRNA %in% rsRNAs_to_annotate) %>%
  select(rsRNA, log2control_median, log2treatment_median) %>%
  pivot_longer(cols = starts_with("log2"), names_to = "Condition", values_to = "Expression") %>%
  mutate(rsRNA = factor(name_mapping[rsRNA], levels = new_names))  # 使用新名称

# 定义颜色从白色到深蓝色
start_color <- rgb(255/255, 255/255, 255/255)  # 白色
end_color <- rgb(0/255, 51/255, 102/255)        # 深蓝色

# 生成渐变
num_colors <- 10
gradient_colors <- colorRampPalette(c(start_color, end_color))(num_colors)

# 绘制热图
heatmap_plot <- ggplot(heatmap_data, aes(x = rsRNA, y = Condition, fill = Expression)) +
  geom_tile(color = "white") +
  scale_fill_gradientn(colors = gradient_colors, 
                       name = "Expression Level") +
  labs(x = "rsRNA", y = "Condition", title = "Heatmap of rsRNA Expression") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 保存热图
ggsave("./bp_rsrna/5bp_heatmap.pdf", plot = heatmap_plot, device = "pdf", dpi = 600, width = 8, height = 6)










# 图中标注名称
miRNAs_to_annotate <- c("")

# Extract the data for miRNAs to annotate
annotations <- merged_data[merged_data$rsRNA %in% miRNAs_to_annotate, ]

# Add annotations for the top 5 miRNAs with the highest p-value (lowest significance)
top_p_values <- merged_data[order(merged_data$lgpadj, decreasing = TRUE), ]
top_p_annotations <- head(top_p_values, 8)

# Combine the annotations: predefined list and the top 5 p-value miRNAs
all_annotations <- rbind(annotations, top_p_annotations)

# Create volcano plot
volcano_plot <- ggplot(merged_data, aes(x = log2fold_change, y = lgpadj)) +
  geom_point(aes(color = ifelse(Q_value < 0.05 & log2fold_change > 0.58, "Up", 
                                ifelse(Q_value < 0.05 & log2fold_change < -0.58, "Down", "Not Sig"))), size = 2) +
  scale_color_manual(values = c("Up" = "#D32F2F", "Down" = "#1976D2", "Not Sig" = "grey70")) +
  labs(x = "Log2 Fold Change", y = "-Log10(padj)", title = "Volcano Plot") +
  theme_minimal(base_size = 14) +
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", size = 0.5) +
  geom_vline(xintercept = c(-0.58, 0.58), linetype = "dashed", color = "black", size = 0.5) +
  geom_point(data = annotations, aes(color = "Highlighted"), size = 3, shape = 1, stroke = 1.5, fill = NA) +  # Highlight miRNAs from predefined list
  geom_text(data = all_annotations, aes(label = rsRNA), vjust = -1, hjust = 0.5, size = 3, color = "black")  # Add text labels for both predefined and top p-value miRNAs

# Display the volcano plot
print(volcano_plot)


ggsave("./bp_rsrna/label_rsrna_volcano_plot3.pdf", plot = volcano_plot, device = "pdf", dpi = 600, width = 16, height = 12)

ggsave("./bp_rsrna/label_rsrna_volcano_plot2.pdf", plot = volcano_plot, device = "pdf", dpi = 600, width = 8, height = 6)
```





## 4. bp_tsrna

2. 标注control_vs_treatment，手动标注
```r
setwd("D:/adult_dep/DE2")
df <- read.csv("./input/rpm_50%_log_tsrna.csv", header = TRUE, row.names = 1)
df <- df[, df[1, ] %in% c(1, 2)]
colnames(df) <- ifelse(df[1, ] == 2, 
                       paste0("control_", colnames(df)), 
                       paste0("treatment_", colnames(df)))
control_columns <- colnames(df)[grepl("^control_", colnames(df))]        

# 删除第一行
df <- df[-1, ]
# 查看修改后的数据框
head(df)
write.csv(df, file = "./input_rename/bp_tsrna.csv")
```

3. 找到每个mirna的treatment_vs_control中位数，计算FC
```bash
cd /mnt/d/adult_dep/DE2
# 记得加“tsrna”
python median.py ./input_rename/bp_tsrna.csv ./bp_tsrna/median.csv
```


② 缺失值替换为NA
```r
setwd("D:/adult_dep/DE2")
data <- read.csv("./input_rename/bp_tsrna.csv") 
data[data == 0] <- NA
write.csv(data, "./bp_tsrna/tsrna_NA.csv", row.names = FALSE)  
```



③ Wilcoxon秩和检验（相关样本）  

秩和检验（Wilcoxon秩和检验）是一种用于比较两组相关或无关样本的非参数检验方法。它不要求数据符合正态分布，因此适用于那些无法满足正态性假设的情况。Wilcoxon秩和检验（相关样本）： 适用于配对设计，比如同一组个体在两个不同时间点的观测。在这种情况下，对每对相关样本进行秩和的比较。    
* 原理：在Wilcoxon秩和检验中，基本思想是将每个组的观测值按大小排序，然后为每个观测值分配一个秩次，最后通过比较秩和来判断两组样本是否来自同一分布。这个检验的零假设是两组样本来自相同的分布。  


```r
# 设置工作目录
setwd("D:/adult_dep/DE2")

# 读取数据
data <- read.csv("./bp_tsrna/tsrna_NA.csv")

# 加载dplyr包
library(dplyr)

# 获取控制组和治疗组的列
control_columns <- grep("^control_", names(data), value = TRUE)
treatment_columns <- grep("^treatment_", names(data), value = TRUE)

# 执行 Wilcoxon 秩和检验
p_values_matrix <- matrix(NA, nrow = length(data$tsRNA), ncol = 3)
colnames(p_values_matrix) <- c("tsRNA", "p_value", "Q_value")

for (i in seq_along(data$tsRNA)) {
  tsRNA <- data$tsRNA[i]
  control_values <- unlist(data[data$tsRNA == tsRNA, control_columns])
  treatment_values <- unlist(data[data$tsRNA == tsRNA, treatment_columns])
  
  # 执行 Wilcoxon 秩和检验
  result <- try(wilcox.test(treatment_values, control_values))

  # 存储 p-value
  if (!inherits(result, "try-error")) {
    p_values_matrix[i, 1] <- tsRNA
    p_values_matrix[i, 2] <- result$p.value
  } else {
    cat("Error for tsRNA:", tsRNA, "\n")
  }
}

# 使用 Benjamini-Hochberg 方法进行多重比较校正
p_values_matrix[, 3] <- p.adjust(p_values_matrix[, 2], method = "BH")

# 输出结果
print(p_values_matrix)
# 将 p_values_matrix 保存到 CSV 文件
write.csv(p_values_matrix, file = "./bp_tsrna/de_values_matrix.csv", row.names = FALSE)
```

④ 合并文件
```r
library(dplyr)

median_data <- read.csv("./bp_tsrna/median.csv")
de_values_data <- read.csv("./bp_tsrna/de_values_matrix.csv")
# 基于 miRNA 列合并数据
merged_data <- inner_join(median_data, de_values_data, by = "tsRNA")

# 查看合并后的数据
print(merged_data)
write.csv(merged_data, file = "./bp_tsrna/result.csv", row.names = FALSE)

# padj 小于 0.05 并且 Log2FC 大于 1（2倍） 或者小于 -1（1/2倍）
diff_gene <- subset(merged_data, Q_value < 0.05 & abs(log2fold_change) > 1)
write.csv(diff_gene, file="./bp_tsrna/sig_q0.05_fc2.csv", quote = F)
```

```r
library(dplyr)
setwd("D:/adult_dep/DE2")

merged_data <- read.csv("./bp_tsrna/result.csv")

diff_gene <- subset(merged_data, Q_value < 0.05 & abs(log2fold_change) > 0.58)
write.csv(diff_gene, file="./bp_tsrna/sig_q0.05_fc15.csv", quote = F)
```



⑤ 绘图
```r
# 加载包
library(ggplot2)
library(pheatmap)

merged_data$lgpadj <- -log10(merged_data$Q_value)

# 绘制火山图
volcano_plot <- ggplot(merged_data, aes(x = log2fold_change, y = lgpadj)) +
  geom_point(aes(color = ifelse(Q_value < 0.05 & log2fold_change > 1, "Up", 
                                ifelse(Q_value < 0.05 & log2fold_change < -1, "Down", "Not Sig"))), size = 2) +
  scale_color_manual(values = c("Up" = "#D32F2F", "Down" = "#1976D2", "Not Sig" = "grey70")) +
  labs(x = "Log2 Fold Change", y = "-Log10(padj)", title = "Volcano Plot") +
  theme_minimal(base_size = 14) +
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(),  # 去掉主网格线
        panel.grid.minor = element_blank(),  # 去掉次网格线
        axis.line = element_line(color = "black")) + # 显示横纵坐标轴
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", size = 0.5) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black", size = 0.5)

ggsave("./bp_tsrna/volcano_plot.pdf", plot = volcano_plot, device = "pdf", dpi = 600, width = 8, height = 6)
```

```r
# 加载包
library(ggplot2)
library(pheatmap)

merged_data <- read.csv("./bp_tsrna/result.csv")
merged_data$lgpadj <- -log10(merged_data$Q_value)

# 绘制火山图
volcano_plot <- ggplot(merged_data, aes(x = log2fold_change, y = lgpadj)) +
  geom_point(aes(color = ifelse(Q_value < 0.05 & log2fold_change > 0.58, "Up", 
                                ifelse(Q_value < 0.05 & log2fold_change < -0.58, "Down", "Not Sig"))), size = 2) +
  scale_color_manual(values = c("Up" = "#D32F2F", "Down" = "#1976D2", "Not Sig" = "grey70")) +
  labs(x = "Log2 Fold Change", y = "-Log10(padj)", title = "Volcano Plot") +
  theme_minimal(base_size = 14) +
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(),  # 去掉主网格线
        panel.grid.minor = element_blank(),  # 去掉次网格线
        axis.line = element_line(color = "black"), # 显示横纵坐标轴
        axis.ticks = element_line(color = "black"),  # 设置坐标轴刻度线颜色
        axis.ticks.length = unit(0.25, "cm")) +  # 设置坐标轴刻度线长度
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", size = 0.5) +
  geom_vline(xintercept = c(-0.58, 0.58), linetype = "dashed", color = "black", size = 0.5)

ggsave("./bp_tsrna/bp_tsrna_volcano_plot.pdf", plot = volcano_plot, device = "pdf", dpi = 600, width = 8, height = 6)
```




* 图中标注

```r
# 图中标注名称
annotations <- data.frame(
  tsRNA = c("ACTTGACCGCTCTGACCA", "ATCCTGCCGACTACGCCA", "GCATTGGTGGTTCAGTGGTAGAATTCT","GCCCGGATAGCTCAGTCGGTAGAG","GCGGGAGACCGGGGTTCGATTCCCCGACGGGGAGC","GTTTCCGTAGTGTAGTGGTTATC"),
  new_label = c("tsRNA-#9", "tsRNA-#17", "tsRNA-#36","tsRNA-#45","tsRNA-#51","tsRNA-#56")
)

# 提取需要标注的 miRNA 数据
annotations_data <- merged_data[merged_data$tsRNA %in% annotations$tsRNA, ]

# 创建火山图
volcano_plot <- ggplot(merged_data, aes(x = log2fold_change, y = lgpadj)) +
  geom_point(aes(color = ifelse(Q_value < 0.05 & log2fold_change > 0.58, "Up", 
                                ifelse(Q_value < 0.05 & log2fold_change < -0.58, "Down", "Not Sig"))), size = 2) +
  scale_color_manual(values = c("Up" = "#D32F2F", "Down" = "#1976D2", "Not Sig" = "grey70")) +
  labs(x = "Log2 Fold Change", y = "-Log10(padj)", title = "Volcano Plot") +
  theme_minimal(base_size = 14) +
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", size = 0.5) +
  geom_vline(xintercept = c(-0.58, 0.58), linetype = "dashed", color = "black", size = 0.5) +
  geom_point(data = annotations_data, aes(color = "Highlighted"), size = 3, shape = 1, stroke = 1.5, fill = NA) +  # 添加圈圈
  geom_text(data = merge(annotations_data, annotations, by = "tsRNA"), aes(label = new_label), 
            vjust = -1, hjust = 0.5, size = 3, color = "black")  # 使用新的标签


ggsave("./bp_tsrna/label_bp_tsrna_volcano_plot.pdf", plot = volcano_plot, device = "pdf", dpi = 600, width = 8, height = 6)








library(ggplot2)
library(pheatmap)

merged_data <- read.csv("./bp_tsrna/result.csv")
merged_data$lgpadj <- -log10(merged_data$Q_value)
# 图中标注名称
miRNAs_to_annotate <- c("")

# Extract the data for miRNAs to annotate
annotations <- merged_data[merged_data$tsRNA %in% miRNAs_to_annotate, ]

# Add annotations for the top 5 miRNAs with the highest p-value (lowest significance)
top_p_values <- merged_data[order(merged_data$lgpadj, decreasing = TRUE), ]
top_p_annotations <- head(top_p_values, 8)

# Combine the annotations: predefined list and the top 5 p-value miRNAs
all_annotations <- rbind(annotations, top_p_annotations)

# Create volcano plot
volcano_plot <- ggplot(merged_data, aes(x = log2fold_change, y = lgpadj)) +
  geom_point(aes(color = ifelse(Q_value < 0.05 & log2fold_change > 0.58, "Up", 
                                ifelse(Q_value < 0.05 & log2fold_change < -0.58, "Down", "Not Sig"))), size = 2) +
  scale_color_manual(values = c("Up" = "#D32F2F", "Down" = "#1976D2", "Not Sig" = "grey70")) +
  labs(x = "Log2 Fold Change", y = "-Log10(padj)", title = "Volcano Plot") +
  theme_minimal(base_size = 14) +
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", size = 0.5) +
  geom_vline(xintercept = c(-0.58, 0.58), linetype = "dashed", color = "black", size = 0.5) +
  geom_point(data = annotations, aes(color = "Highlighted"), size = 3, shape = 1, stroke = 1.5, fill = NA) +  # Highlight miRNAs from predefined list
  geom_text(data = all_annotations, aes(label = tsRNA), vjust = -1, hjust = 0.5, size = 3, color = "black")  # Add text labels for both predefined and top p-value miRNAs

# Display the volcano plot
print(volcano_plot)


ggsave("./bp_tsrna/label_tsrna_volcano_plot3.pdf", plot = volcano_plot, device = "pdf", dpi = 600, width = 16, height = 12)

ggsave("./bp_tsrna/label_tsrna_volcano_plot2.pdf", plot = volcano_plot, device = "pdf", dpi = 600, width = 8, height = 6)
```

```r
# 绘制热图
library(ggplot2)
library(tidyr)
library(dplyr)

# 读取数据
merged_data <- read.csv("./bp_tsrna/result.csv")

# 指定需要标注的 tsRNA 和对应的新名称
tsRNAs_to_annotate <- c("ACTTGACCGCTCTGACCA", "ATCCTGCCGACTACGCCA", "GCATTGGTGGTTCAGTGGTAGAATTCT","GCCCGGATAGCTCAGTCGGTAGAG","GCGGGAGACCGGGGTTCGATTCCCCGACGGGGAGC","GTTTCCGTAGTGTAGTGGTTATC")

new_names <- c("tsRNA-#9", "tsRNA-#17", "tsRNA-#36","tsRNA-#45","tsRNA-#51","tsRNA-#56")

# 创建一个命名映射
name_mapping <- setNames(new_names, tsRNAs_to_annotate)

# 提取需要的 tsRNA 数据并设置顺序
heatmap_data <- merged_data %>%
  filter(tsRNA %in% tsRNAs_to_annotate) %>%
  select(tsRNA, log2control_median, log2treatment_median) %>%
  pivot_longer(cols = starts_with("log2"), names_to = "Condition", values_to = "Expression") %>%
  mutate(tsRNA = factor(name_mapping[tsRNA], levels = new_names))  # 使用新名称

# 定义颜色从白色到深蓝色
start_color <- rgb(255/255, 255/255, 255/255)  # 白色
end_color <- rgb(0/255, 51/255, 102/255)        # 深蓝色

# 生成渐变
num_colors <- 10
gradient_colors <- colorRampPalette(c(start_color, end_color))(num_colors)

# 绘制热图
heatmap_plot <- ggplot(heatmap_data, aes(x = tsRNA, y = Condition, fill = Expression)) +
  geom_tile(color = "white") +
  scale_fill_gradientn(colors = gradient_colors, 
                       name = "Expression Level") +
  labs(x = "tsRNA", y = "Condition", title = "Heatmap of tsRNA Expression") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 保存热图
ggsave("./bp_tsrna/5tsRNA_heatmap.pdf", plot = heatmap_plot, device = "pdf", dpi = 600, width = 8, height = 6)
```



























# hc+36other_MDD

## 1. hc+36other_mirna

2. 标注control_vs_treatment，手动标注
```r
setwd("D:/adult_dep/DE2")
df <- read.csv("./input/rpm_50%_log_mirna.csv", header = TRUE, row.names = 1)
df <- df[, df[1, ] %in% c(1, 0, 3)]
colnames(df) <- ifelse(df[1, ] %in% c(0, 3), 
                       paste0("control_", colnames(df)), 
                       paste0("treatment_", colnames(df)))
control_columns <- colnames(df)[grepl("^control_", colnames(df))] 
treatment_columns <- colnames(df)[grepl("^treatment_", colnames(df))]        


# 删除第一行
df <- df[-1, ]
# 查看修改后的数据框
head(df)
write.csv(df, file = "./input_rename/hc+36other_mirna.csv")
```

3. 找到每个mirna的treatment_vs_control中位数，计算FC
```bash
cd /mnt/d/adult_dep/DE2
# 记得加“mirna”
python median_mirna.py ./input_rename/hc+36other_mirna.csv ./hc+36other_mirna/median.csv
```


② 缺失值替换为NA
```r
setwd("D:/adult_dep/DE2")
data <- read.csv("./input_rename/hc+36other_mirna.csv") 
data[data == 0] <- NA
write.csv(data, "./hc+36other_mirna/mirna_NA.csv", row.names = FALSE)  
```



③ Wilcoxon秩和检验（相关样本）  

秩和检验（Wilcoxon秩和检验）是一种用于比较两组相关或无关样本的非参数检验方法。它不要求数据符合正态分布，因此适用于那些无法满足正态性假设的情况。Wilcoxon秩和检验（相关样本）： 适用于配对设计，比如同一组个体在两个不同时间点的观测。在这种情况下，对每对相关样本进行秩和的比较。    
* 原理：在Wilcoxon秩和检验中，基本思想是将每个组的观测值按大小排序，然后为每个观测值分配一个秩次，最后通过比较秩和来判断两组样本是否来自同一分布。这个检验的零假设是两组样本来自相同的分布。  


```r
# 设置工作目录
setwd("D:/adult_dep/DE2")

# 读取数据
data <- read.csv("./hc+36other_mirna/mirna_NA.csv")

# 加载dplyr包
library(dplyr)

# 获取控制组和治疗组的列
control_columns <- grep("^control_", names(data), value = TRUE)
treatment_columns <- grep("^treatment_", names(data), value = TRUE)


# 执行 Wilcoxon 秩和检验
p_values_matrix <- matrix(NA, nrow = length(data$miRNA), ncol = 3)
colnames(p_values_matrix) <- c("miRNA", "p_value", "Q_value")

for (i in seq_along(data$miRNA)) {
  miRNA <- data$miRNA[i]
  control_values <- unlist(data[data$miRNA == miRNA, control_columns])
  treatment_values <- unlist(data[data$miRNA == miRNA, treatment_columns])
  
  # 执行 Wilcoxon 秩和检验
  result <- try(wilcox.test(treatment_values, control_values))

  # 存储 p-value
  if (!inherits(result, "try-error")) {
    p_values_matrix[i, 1] <- miRNA
    p_values_matrix[i, 2] <- result$p.value
  } else {
    cat("Error for miRNA:", miRNA, "\n")
  }
}

# 使用 Benjamini-Hochberg 方法进行多重比较校正
p_values_matrix[, 3] <- p.adjust(p_values_matrix[, 2], method = "BH")

# 输出结果
print(p_values_matrix)
# 将 p_values_matrix 保存到 CSV 文件
write.csv(p_values_matrix, file = "./hc+36other_mirna/de_values_matrix.csv", row.names = FALSE)
```

④ 合并文件
```r
library(dplyr)

median_data <- read.csv("./hc+36other_mirna/median.csv")
de_values_data <- read.csv("./hc+36other_mirna/de_values_matrix.csv")
# 基于 miRNA 列合并数据
merged_data <- inner_join(median_data, de_values_data, by = "miRNA")

# 查看合并后的数据
print(merged_data)
write.csv(merged_data, file = "./hc+36other_mirna/result.csv", row.names = FALSE)

diff_gene <- subset(merged_data, Q_value < 0.05 & abs(log2fold_change) > 0.58)
write.csv(diff_gene, file="./hc+36other_mirna/sig_q0.05_fc15.csv", quote = F)
```



⑤ 绘图
```r
# 加载包
library(ggplot2)
library(pheatmap)
setwd("D:/adult_dep/DE2")

merged_data <- read.csv("./hc+36other_mirna/result.csv")
merged_data$lgpadj <- -log10(merged_data$Q_value)

# Define categories for points
merged_data$category <- ifelse(merged_data$Q_value < 0.05 & merged_data$log2fold_change > 0.58, "Up", 
                                ifelse(merged_data$Q_value < 0.05 & merged_data$log2fold_change < -0.58, "Down", "Not Sig"))

# Calculate counts for each category
category_counts <- table(merged_data$category)

# Create a volcano plot
volcano_plot <- ggplot(merged_data, aes(x = log2fold_change, y = lgpadj)) +
  geom_point(aes(color = category), size = 2) +
  scale_color_manual(values = c("Up" = "#D32F2F", "Down" = "#1976D2", "Not Sig" = "grey70")) +
  labs(x = "Log2 Fold Change", y = "-Log10(padj)", title = "Volcano Plot") +
  theme_minimal(base_size = 14) +
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        axis.line = element_line(color = "black"),  # Show axis lines
        axis.ticks = element_line(color = "black"),  # Show axis ticks
        axis.ticks.length = unit(0.25, "cm")) +  # Set tick length
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", size = 0.5) +
  geom_vline(xintercept = c(-0.58, 0.58), linetype = "dashed", color = "black", size = 0.5) +
  
  # Add category count labels
  geom_text(aes(x = 2, y = 30, label = paste("Up:", category_counts["Up"])), color = "#D32F2F", size = 5, hjust = 0) +
  geom_text(aes(x = 2, y = 28, label = paste("Down:", category_counts["Down"])), color = "#1976D2", size = 5, hjust = 0) +
  geom_text(aes(x = 2, y = 26, label = paste("Not Sig:", category_counts["Not Sig"])), color = "grey70", size = 5, hjust = 0)

# Save the plot
ggsave("./hc+36other_mirna/hc+36other_mirna_volcano_plot.pdf", plot = volcano_plot, device = "pdf", dpi = 600, width = 8, height = 6)
```

```r
# 图中标注名称
# miRNAs to annotate (predefined list)
miRNAs_to_annotate <- c("hsa-mir-766_hsa-miR-766-3p", "hsa-mir-1277_hsa-miR-1277-5p", 
                         "hsa-mir-548k_hsa-miR-548k", "hsa-mir-151b_hsa-miR-151b", 
                         "hsa-mir-4433b_hsa-miR-4433b-3p", "hsa-mir-584_hsa-miR-584-5p")

# Extract the data for miRNAs to annotate
annotations <- merged_data[merged_data$miRNA %in% miRNAs_to_annotate, ]

# Add annotations for the top 5 miRNAs with the highest p-value (lowest significance)
top_p_values <- merged_data[order(merged_data$lgpadj, decreasing = TRUE), ]
top_p_annotations <- head(top_p_values, 5)

# Combine the annotations: predefined list and the top 5 p-value miRNAs
all_annotations <- rbind(annotations, top_p_annotations)

# Create volcano plot
volcano_plot <- ggplot(merged_data, aes(x = log2fold_change, y = lgpadj)) +
  geom_point(aes(color = ifelse(Q_value < 0.05 & log2fold_change > 0.58, "Up", 
                                ifelse(Q_value < 0.05 & log2fold_change < -0.58, "Down", "Not Sig"))), size = 2) +
  scale_color_manual(values = c("Up" = "#D32F2F", "Down" = "#1976D2", "Not Sig" = "grey70")) +
  labs(x = "Log2 Fold Change", y = "-Log10(padj)", title = "Volcano Plot") +
  theme_minimal(base_size = 14) +
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", size = 0.5) +
  geom_vline(xintercept = c(-0.58, 0.58), linetype = "dashed", color = "black", size = 0.5) +
  geom_point(data = annotations, aes(color = "Highlighted"), size = 3, shape = 1, stroke = 1.5, fill = NA) +  # Highlight miRNAs from predefined list
  geom_text(data = all_annotations, aes(label = miRNA), vjust = -1, hjust = 0.5, size = 3, color = "black")  # Add text labels for both predefined and top p-value miRNAs

# Display the volcano plot
print(volcano_plot)


ggsave("./hc+36other_mirna/label_hc+36other_mirna_volcano_plot3.pdf", plot = volcano_plot, device = "pdf", dpi = 600, width = 16, height = 12)


ggsave("./hc+36other_mirna/label_hc+36other_mirna_volcano_plot2.pdf", plot = volcano_plot, device = "pdf", dpi = 600, width = 8, height = 6)

```

```r
# 绘制热图
library(ggplot2)
library(tidyr)
library(dplyr)

# 读取数据
merged_data <- read.csv("./other_mirna/result.csv")

# 指定要绘制热图的 miRNA
miRNAs_to_annotate <- c("hsa-mir-101-1_hsa-miR-101-3p", "hsa-mir-140_hsa-miR-140-3p", "hsa-mir-25_hsa-miR-25-3p",
                        "hsa-mir-20a_hsa-miR-20a-5p","hsa-mir-7-3_hsa-miR-7-5p","hsa-mir-15b_hsa-miR-15b-3p", "hsa-mir-19b-1_hsa-miR-19b-3p")

# 提取需要的 miRNA 数据并设置顺序
heatmap_data <- merged_data %>%
  filter(miRNA %in% miRNAs_to_annotate) %>%
  select(miRNA, log2control_median, log2treatment_median) %>%
  pivot_longer(cols = starts_with("log2"), names_to = "Condition", values_to = "Expression") %>%
  mutate(miRNA = factor(miRNA, levels = miRNAs_to_annotate))  # 设置因子顺序

# 定义颜色从白色到深蓝色
start_color <- rgb(255/255, 255/255, 255/255)  # 白色
end_color <- rgb(0/255, 51/255, 102/255)        # 深蓝色

# 生成渐变
num_colors <- 10
gradient_colors <- colorRampPalette(c(start_color, end_color))(num_colors)

# 绘制热图
heatmap_plot <- ggplot(heatmap_data, aes(x = miRNA, y = Condition, fill = Expression)) +
  geom_tile(color = "white") +
  scale_fill_gradientn(colors = gradient_colors, 
                       name = "Expression Level") +
  labs(x = "miRNA", y = "Condition", title = "Heatmap of miRNA Expression") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 保存热图
ggsave("./other_mirna/5mirna_heatmap.pdf", plot = heatmap_plot, device = "pdf", dpi = 600, width = 8, height = 6)
```




## 2. hc+36other_pirna

2. 标注control_vs_treatment，手动标注
```r
setwd("D:/adult_dep/DE2")
df <- read.csv("./input/rpm_50%_log_pirna.csv", header = TRUE, row.names = 1)
df <- df[, df[1, ] %in% c(1, 0, 3)]
colnames(df) <- ifelse(df[1, ] %in% c(0, 3), 
                       paste0("control_", colnames(df)), 
                       paste0("treatment_", colnames(df)))
control_columns <- colnames(df)[grepl("^control_", colnames(df))]        

# 删除第一行
df <- df[-1, ]
# 查看修改后的数据框
head(df)
write.csv(df, file = "./input_rename/hc+36other_pirna.csv")
```

3. 找到每个mirna的treatment_vs_control中位数，计算FC
```bash
cd /mnt/d/adult_dep/DE2
# 记得加“pirna”
python median_pirna.py ./input_rename/hc+36other_pirna.csv ./hc+36other_pirna/median.csv
```


② 缺失值替换为NA
```r
setwd("D:/adult_dep/DE2")
data <- read.csv("./input_rename/hc+36other_pirna.csv") 
data[data == 0] <- NA
write.csv(data, "./hc+36other_pirna/pirna_NA.csv", row.names = FALSE)  
```



③ Wilcoxon秩和检验（相关样本）  

秩和检验（Wilcoxon秩和检验）是一种用于比较两组相关或无关样本的非参数检验方法。它不要求数据符合正态分布，因此适用于那些无法满足正态性假设的情况。Wilcoxon秩和检验（相关样本）： 适用于配对设计，比如同一组个体在两个不同时间点的观测。在这种情况下，对每对相关样本进行秩和的比较。    
* 原理：在Wilcoxon秩和检验中，基本思想是将每个组的观测值按大小排序，然后为每个观测值分配一个秩次，最后通过比较秩和来判断两组样本是否来自同一分布。这个检验的零假设是两组样本来自相同的分布。  


```r
# 设置工作目录
setwd("D:/adult_dep/DE2")

# 读取数据
data <- read.csv("./hc+36other_pirna/pirna_NA.csv")

# 加载dplyr包
library(dplyr)

# 获取控制组和治疗组的列
control_columns <- grep("^control_", names(data), value = TRUE)
treatment_columns <- grep("^treatment_", names(data), value = TRUE)
# 835

# 执行 Wilcoxon 秩和检验
p_values_matrix <- matrix(NA, nrow = length(data$piRNA), ncol = 3)
colnames(p_values_matrix) <- c("piRNA", "p_value", "Q_value")

for (i in seq_along(data$piRNA)) {
  piRNA <- data$piRNA[i]
  control_values <- unlist(data[data$piRNA == piRNA, control_columns])
  treatment_values <- unlist(data[data$piRNA == piRNA, treatment_columns])
  
  # 执行 Wilcoxon 秩和检验
  result <- try(wilcox.test(treatment_values, control_values))

  # 存储 p-value
  if (!inherits(result, "try-error")) {
    p_values_matrix[i, 1] <- piRNA
    p_values_matrix[i, 2] <- result$p.value
  } else {
    cat("Error for piRNA:", piRNA, "\n")
  }
}

# 使用 Benjamini-Hochberg 方法进行多重比较校正
p_values_matrix[, 3] <- p.adjust(p_values_matrix[, 2], method = "BH")

# 输出结果
print(p_values_matrix)
# 将 p_values_matrix 保存到 CSV 文件
write.csv(p_values_matrix, file = "./hc+36other_pirna/de_values_matrix.csv", row.names = FALSE)
```

④ 合并文件
```r
library(dplyr)

median_data <- read.csv("./hc+36other_pirna/median.csv")
de_values_data <- read.csv("./hc+36other_pirna/de_values_matrix.csv")
# 基于 miRNA 列合并数据
merged_data <- inner_join(median_data, de_values_data, by = "piRNA")
# 查看合并后的数据
print(merged_data)
write.csv(merged_data, file = "./hc+36other_pirna/result.csv", row.names = FALSE)

diff_gene <- subset(merged_data, Q_value < 0.05 & abs(log2fold_change) > 0.58)
write.csv(diff_gene, file="./hc+36other_pirna/sig_q0.05_fc15.csv", quote = F)
```


⑤ 绘图
```r
# 加载包
library(ggplot2)
library(pheatmap)

merged_data$lgpadj <- -log10(merged_data$Q_value)

# 绘制火山图
volcano_plot <- ggplot(merged_data, aes(x = log2fold_change, y = lgpadj)) +
  geom_point(aes(color = ifelse(Q_value < 0.05 & log2fold_change > 1, "Up", 
                                ifelse(Q_value < 0.05 & log2fold_change < -1, "Down", "Not Sig"))), size = 2) +
  scale_color_manual(values = c("Up" = "#D32F2F", "Down" = "#1976D2", "Not Sig" = "grey70")) +
  labs(x = "Log2 Fold Change", y = "-Log10(padj)", title = "Volcano Plot") +
  theme_minimal(base_size = 14) +
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(),  # 去掉主网格线
        panel.grid.minor = element_blank(),  # 去掉次网格线
        axis.line = element_line(color = "black")) + # 显示横纵坐标轴
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", size = 0.5) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black", size = 0.5)

ggsave("./other_pirna/volcano_plot.pdf", plot = volcano_plot, device = "pdf", dpi = 600, width = 8, height = 6)
```

```r
# 加载包
library(ggplot2)
library(pheatmap)

merged_data <- read.csv("./hc+36other_pirna/result.csv")
merged_data$lgpadj <- -log10(merged_data$Q_value)

# Define categories for points
merged_data$category <- ifelse(merged_data$Q_value < 0.05 & merged_data$log2fold_change > 0.58, "Up", 
                                ifelse(merged_data$Q_value < 0.05 & merged_data$log2fold_change < -0.58, "Down", "Not Sig"))

# Calculate counts for each category
category_counts <- table(merged_data$category)

# Create a volcano plot
volcano_plot <- ggplot(merged_data, aes(x = log2fold_change, y = lgpadj)) +
  geom_point(aes(color = category), size = 2) +
  scale_color_manual(values = c("Up" = "#D32F2F", "Down" = "#1976D2", "Not Sig" = "grey70")) +
  labs(x = "Log2 Fold Change", y = "-Log10(padj)", title = "Volcano Plot") +
  theme_minimal(base_size = 14) +
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        axis.line = element_line(color = "black"),  # Show axis lines
        axis.ticks = element_line(color = "black"),  # Show axis ticks
        axis.ticks.length = unit(0.25, "cm")) +  # Set tick length
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", size = 0.5) +
  geom_vline(xintercept = c(-0.58, 0.58), linetype = "dashed", color = "black", size = 0.5) +
  
  # Add category count labels
  geom_text(aes(x = 1, y = 30, label = paste("Up:", category_counts["Up"])), color = "#D32F2F", size = 5, hjust = 0) +
  geom_text(aes(x = 1, y = 28, label = paste("Down:", category_counts["Down"])), color = "#1976D2", size = 5, hjust = 0) +
  geom_text(aes(x = 1, y = 26, label = paste("Not Sig:", category_counts["Not Sig"])), color = "grey70", size = 5, hjust = 0)

# Save the plot
ggsave("./hc+36other_pirna/hc+36other_pirna_volcano_plot.pdf", plot = volcano_plot, device = "pdf", dpi = 600, width = 8, height = 6)
```
* 图中标注

```r
# 图中标注名称
miRNAs_to_annotate <- c("piR-hsa-187496", "piR-hsa-369497","piR-hsa-188007","piR-hsa-4419732","piR-hsa-334268")

# Extract the data for miRNAs to annotate
annotations <- merged_data[merged_data$piRNA %in% miRNAs_to_annotate, ]

# Add annotations for the top 5 miRNAs with the highest p-value (lowest significance)
top_p_values <- merged_data[order(merged_data$lgpadj, decreasing = TRUE), ]
top_p_annotations <- head(top_p_values, 3)

# Combine the annotations: predefined list and the top 5 p-value miRNAs
all_annotations <- rbind(annotations, top_p_annotations)

# Create volcano plot
volcano_plot <- ggplot(merged_data, aes(x = log2fold_change, y = lgpadj)) +
  geom_point(aes(color = ifelse(Q_value < 0.05 & log2fold_change > 0.58, "Up", 
                                ifelse(Q_value < 0.05 & log2fold_change < -0.58, "Down", "Not Sig"))), size = 2) +
  scale_color_manual(values = c("Up" = "#D32F2F", "Down" = "#1976D2", "Not Sig" = "grey70")) +
  labs(x = "Log2 Fold Change", y = "-Log10(padj)", title = "Volcano Plot") +
  theme_minimal(base_size = 14) +
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", size = 0.5) +
  geom_vline(xintercept = c(-0.58, 0.58), linetype = "dashed", color = "black", size = 0.5) +
  geom_point(data = annotations, aes(color = "Highlighted"), size = 3, shape = 1, stroke = 1.5, fill = NA) +  # Highlight miRNAs from predefined list
  geom_text(data = all_annotations, aes(label = piRNA), vjust = -1, hjust = 0.5, size = 3, color = "black")  # Add text labels for both predefined and top p-value miRNAs

# Display the volcano plot
print(volcano_plot)


ggsave("./hc+36other_pirna/label_hc+36other_pirna_volcano_plot3.pdf", plot = volcano_plot, device = "pdf", dpi = 600, width = 16, height = 12)

ggsave("./hc+36other_pirna/label_hc+36other_pirna_volcano_plot2.pdf", plot = volcano_plot, device = "pdf", dpi = 600, width = 8, height = 6)

```

```r
# 绘制热图
library(ggplot2)
library(tidyr)
library(dplyr)

# 读取数据
merged_data <- read.csv("./other_pirna/result.csv")

# 指定要绘制热图的 miRNA
piRNAs_to_annotate <- c("piR-hsa-168711", "piR-hsa-2482471","piR-hsa-334268","piR-hsa-35389")

# 提取需要的 piRNA 数据并设置顺序
heatmap_data <- merged_data %>%
  filter(piRNA %in% piRNAs_to_annotate) %>%
  select(piRNA, log2control_median, log2treatment_median) %>%
  pivot_longer(cols = starts_with("log2"), names_to = "Condition", values_to = "Expression") %>%
  mutate(piRNA = factor(piRNA, levels = piRNAs_to_annotate))  # 设置因子顺序

# 定义颜色从白色到深蓝色
start_color <- rgb(255/255, 255/255, 255/255)  # 白色
end_color <- rgb(0/255, 51/255, 102/255)        # 深蓝色

# 生成渐变
num_colors <- 10
gradient_colors <- colorRampPalette(c(start_color, end_color))(num_colors)

# 绘制热图
heatmap_plot <- ggplot(heatmap_data, aes(x = piRNA, y = Condition, fill = Expression)) +
  geom_tile(color = "white") +
  scale_fill_gradientn(colors = gradient_colors, 
                       name = "Expression Level") +
  labs(x = "piRNA", y = "Condition", title = "Heatmap of piRNA Expression") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 保存热图
ggsave("./other_pirna/5piRNA_heatmap.pdf", plot = heatmap_plot, device = "pdf", dpi = 600, width = 8, height = 6)
```




## 3. hc+36other_rsrna

2. 标注control_vs_treatment，手动标注
```r
setwd("D:/adult_dep/DE2")
df <- read.csv("./input/rpm_50%_log_rsrna.csv", header = TRUE, row.names = 1)
df <- df[, df[1, ] %in% c(1, 0, 3)]
colnames(df) <- ifelse(df[1, ] %in% c(0, 3), 
                       paste0("control_", colnames(df)), 
                       paste0("treatment_", colnames(df)))
control_columns <- colnames(df)[grepl("^control_", colnames(df))]        
treatment_columns <- colnames(df)[grepl("^treatment_", colnames(df))]        

# 删除第一行
df <- df[-1, ]
# 查看修改后的数据框
head(df)
write.csv(df, file = "./input_rename/hc+36other_rsrna.csv")
```

3. 找到每个mirna的treatment_vs_control中位数，计算FC
```bash
cd /mnt/d/adult_dep/DE2
# 记得加“rsrna”
python median_rsrna.py ./input_rename/hc+36other_rsrna.csv ./hc+36other_rsrna/median.csv
```


② 缺失值替换为NA
```r
setwd("D:/adult_dep/DE2")
data <- read.csv("./input_rename/hc+36other_rsrna.csv") 
data[data == 0] <- NA
write.csv(data, "./hc+36other_rsrna/rsrna_NA.csv", row.names = FALSE)  
```



③ Wilcoxon秩和检验（相关样本）  

秩和检验（Wilcoxon秩和检验）是一种用于比较两组相关或无关样本的非参数检验方法。它不要求数据符合正态分布，因此适用于那些无法满足正态性假设的情况。Wilcoxon秩和检验（相关样本）： 适用于配对设计，比如同一组个体在两个不同时间点的观测。在这种情况下，对每对相关样本进行秩和的比较。    
* 原理：在Wilcoxon秩和检验中，基本思想是将每个组的观测值按大小排序，然后为每个观测值分配一个秩次，最后通过比较秩和来判断两组样本是否来自同一分布。这个检验的零假设是两组样本来自相同的分布。  


```r
# 设置工作目录
setwd("D:/adult_dep/DE2")

# 读取数据
data <- read.csv("./hc+36other_rsrna/rsrna_NA.csv")

# 加载dplyr包
library(dplyr)

# 获取控制组和治疗组的列
control_columns <- grep("^control_", names(data), value = TRUE)
treatment_columns <- grep("^treatment_", names(data), value = TRUE)

# 执行 Wilcoxon 秩和检验
p_values_matrix <- matrix(NA, nrow = length(data$rsRNA), ncol = 3)
colnames(p_values_matrix) <- c("rsRNA", "p_value", "Q_value")

for (i in seq_along(data$rsRNA)) {
  rsRNA <- data$rsRNA[i]
  control_values <- unlist(data[data$rsRNA == rsRNA, control_columns])
  treatment_values <- unlist(data[data$rsRNA == rsRNA, treatment_columns])
  
  # 执行 Wilcoxon 秩和检验
  result <- try(wilcox.test(treatment_values, control_values))

  # 存储 p-value
  if (!inherits(result, "try-error")) {
    p_values_matrix[i, 1] <- rsRNA
    p_values_matrix[i, 2] <- result$p.value
  } else {
    cat("Error for rsRNA:", rsRNA, "\n")
  }
}

# 使用 Benjamini-Hochberg 方法进行多重比较校正
p_values_matrix[, 3] <- p.adjust(p_values_matrix[, 2], method = "BH")

# 输出结果
print(p_values_matrix)
write.csv(p_values_matrix, file = "./hc+36other_rsrna/de_values_matrix.csv", row.names = FALSE)
```

④ 合并文件
```r
library(dplyr)

median_data <- read.csv("./hc+36other_rsrna/median.csv")
de_values_data <- read.csv("./hc+36other_rsrna/de_values_matrix.csv")
# 基于 miRNA 列合并数据
merged_data <- inner_join(median_data, de_values_data, by = "rsRNA")

# 查看合并后的数据
print(merged_data)
write.csv(merged_data, file = "./hc+36other_rsrna/result.csv", row.names = FALSE)

diff_gene <- subset(merged_data, Q_value < 0.05 & abs(log2fold_change) > 0.58)
write.csv(diff_gene, file="./hc+36other_rsrna/sig_q0.05_fc15.csv", quote = F)
```

⑤ 绘图
```r
# 加载包
library(ggplot2)
library(pheatmap)

merged_data$lgpadj <- -log10(merged_data$Q_value)

# 绘制火山图
volcano_plot <- ggplot(merged_data, aes(x = log2fold_change, y = lgpadj)) +
  geom_point(aes(color = ifelse(Q_value < 0.05 & log2fold_change > 1, "Up", 
                                ifelse(Q_value < 0.05 & log2fold_change < -1, "Down", "Not Sig"))), size = 2) +
  scale_color_manual(values = c("Up" = "#D32F2F", "Down" = "#1976D2", "Not Sig" = "grey70")) +
  labs(x = "Log2 Fold Change", y = "-Log10(padj)", title = "Volcano Plot") +
  theme_minimal(base_size = 14) +
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(),  # 去掉主网格线
        panel.grid.minor = element_blank(),  # 去掉次网格线
        axis.line = element_line(color = "black")) + # 显示横纵坐标轴
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", size = 0.5) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black", size = 0.5)

ggsave("./other_rsrna/volcano_plot.pdf", plot = volcano_plot, device = "pdf", dpi = 600, width = 8, height = 6)
```
```r
# 加载包
library(ggplot2)
library(pheatmap)

merged_data <- read.csv("./hc+36other_rsrna/result.csv")
merged_data$lgpadj <- -log10(merged_data$Q_value)

# Define categories for points
merged_data$category <- ifelse(merged_data$Q_value < 0.05 & merged_data$log2fold_change > 0.58, "Up", 
                                ifelse(merged_data$Q_value < 0.05 & merged_data$log2fold_change < -0.58, "Down", "Not Sig"))

# Calculate counts for each category
category_counts <- table(merged_data$category)

# Create a volcano plot
volcano_plot <- ggplot(merged_data, aes(x = log2fold_change, y = lgpadj)) +
  geom_point(aes(color = category), size = 2) +
  scale_color_manual(values = c("Up" = "#D32F2F", "Down" = "#1976D2", "Not Sig" = "grey70")) +
  labs(x = "Log2 Fold Change", y = "-Log10(padj)", title = "Volcano Plot") +
  theme_minimal(base_size = 14) +
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        axis.line = element_line(color = "black"),  # Show axis lines
        axis.ticks = element_line(color = "black"),  # Show axis ticks
        axis.ticks.length = unit(0.25, "cm")) +  # Set tick length
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", size = 0.5) +
  geom_vline(xintercept = c(-0.58, 0.58), linetype = "dashed", color = "black", size = 0.5) +
  
  # Add category count labels
  geom_text(aes(x = 1, y = 30, label = paste("Up:", category_counts["Up"])), color = "#D32F2F", size = 5, hjust = 0) +
  geom_text(aes(x = 1, y = 28, label = paste("Down:", category_counts["Down"])), color = "#1976D2", size = 5, hjust = 0) +
  geom_text(aes(x = 1, y = 26, label = paste("Not Sig:", category_counts["Not Sig"])), color = "grey70", size = 5, hjust = 0)

# Save the plot
ggsave("./hc+36other_rsrna/hc+36other_rsrna_volcano_plot.pdf", plot = volcano_plot, device = "pdf", dpi = 600, width = 8, height = 6)
```
* 图中标注

```r
# 图中标注名称
annotations <- data.frame(
  rsRNA = c("ACGGACCAAGGAGTCTAACACGTGCGCG", "GGAGTCTAACACGTGCGCG"),
  new_label = c("rsRNA-#29", "rsRNA-#182")
)

# 提取需要标注的 miRNA 数据
annotations_data <- merged_data[merged_data$rsRNA %in% annotations$rsRNA, ]

# 创建火山图
volcano_plot <- ggplot(merged_data, aes(x = log2fold_change, y = lgpadj)) +
  geom_point(aes(color = ifelse(Q_value < 0.05 & log2fold_change > 0.58, "Up", 
                                ifelse(Q_value < 0.05 & log2fold_change < -0.58, "Down", "Not Sig"))), size = 2) +
  scale_color_manual(values = c("Up" = "#D32F2F", "Down" = "#1976D2", "Not Sig" = "grey70")) +
  labs(x = "Log2 Fold Change", y = "-Log10(padj)", title = "Volcano Plot") +
  theme_minimal(base_size = 14) +
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", size = 0.5) +
  geom_vline(xintercept = c(-0.58, 0.58), linetype = "dashed", color = "black", size = 0.5) +
  geom_point(data = annotations_data, aes(color = "Highlighted"), size = 3, shape = 1, stroke = 1.5, fill = NA) +  # 添加圈圈
  geom_text(data = merge(annotations_data, annotations, by = "rsRNA"), aes(label = new_label), 
            vjust = -1, hjust = 0.5, size = 3, color = "black")  # 使用新的标签


ggsave("./other_rsrna/label_other_rsrna_volcano_plot.pdf", plot = volcano_plot, device = "pdf", dpi = 600, width = 8, height = 6)




# 图中标注名称
annotations <- data.frame(
  rsRNA = c("CTGGCACGGTGAAGAGACATGAGAGGTGTAGAATA", "ACCAAGGAGTCTAACACGTGCGCG", "GACCAAGGAGTCTAACACGTGCGCG", "CGGACCAAGGAGTCTAACACGTGCGCG", "CAAGGAGTCTAACACGTGCGCG", "TGGTGGTGCATGGCCGTTCTT", "GTGGTGCATGGCCGTTCTT"),
  new_label = c("rsRNA-#134", "rsRNA-#22", "rsRNA-#146", "rsRNA-#122", "rsRNA-#69", "rsRNA-#267", "rsRNA-#216")
)

# 提取需要标注的 miRNA 数据
annotations <- merged_data[merged_data$rsRNA %in% annotations$rsRNA, ]

# Add annotations for the top 5 miRNAs with the highest p-value (lowest significance)
top_p_values <- merged_data[order(merged_data$lgpadj, decreasing = TRUE), ]
top_p_annotations <- head(top_p_values, 3)

# Combine the annotations: predefined list and the top 5 p-value miRNAs
all_annotations <- rbind(annotations, top_p_annotations)

# Create volcano plot
volcano_plot <- ggplot(merged_data, aes(x = log2fold_change, y = lgpadj)) +
  geom_point(aes(color = ifelse(Q_value < 0.05 & log2fold_change > 0.58, "Up", 
                                ifelse(Q_value < 0.05 & log2fold_change < -0.58, "Down", "Not Sig"))), size = 2) +
  scale_color_manual(values = c("Up" = "#D32F2F", "Down" = "#1976D2", "Not Sig" = "grey70")) +
  labs(x = "Log2 Fold Change", y = "-Log10(padj)", title = "Volcano Plot") +
  theme_minimal(base_size = 14) +
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", size = 0.5) +
  geom_vline(xintercept = c(-0.58, 0.58), linetype = "dashed", color = "black", size = 0.5) +
  geom_point(data = annotations, aes(color = "Highlighted"), size = 3, shape = 1, stroke = 1.5, fill = NA) +  # Highlight miRNAs from predefined list
  geom_text(data = all_annotations, aes(label = rsRNA), vjust = -1, hjust = 0.5, size = 3, color = "black")  # Add text labels for both predefined and top p-value miRNAs

# Display the volcano plot
print(volcano_plot)


ggsave("./hc+36other_rsrna/label_hc+36other_rsrna_volcano_plot3.pdf", plot = volcano_plot, device = "pdf", dpi = 600, width = 16, height = 12)

ggsave("./hc+36other_rsrna/label_hc+36other_rsrna_volcano_plot2.pdf", plot = volcano_plot, device = "pdf", dpi = 600, width = 8, height = 6)
```














## 4. hc+36other_tsrna

2. 标注control_vs_treatment，手动标注
```r
setwd("D:/adult_dep/DE2")
df <- read.csv("./input/rpm_50%_log_tsrna.csv", header = TRUE, row.names = 1)
df <- df[, df[1, ] %in% c(1, 0, 3)]
colnames(df) <- ifelse(df[1, ] %in% c(0, 3), 
                       paste0("control_", colnames(df)), 
                       paste0("treatment_", colnames(df)))
control_columns <- colnames(df)[grepl("^control_", colnames(df))] 
treatment_columns <- colnames(df)[grepl("^treatment_", colnames(df))]        


# 删除第一行
df <- df[-1, ]
# 查看修改后的数据框
head(df)
write.csv(df, file = "./input_rename/hc+36other_tsrna.csv")
```

3. 找到每个mirna的treatment_vs_control中位数，计算FC
```bash
cd /mnt/d/adult_dep/DE2
# 记得加“tsRNA”
python median_tsrna.py ./input_rename/hc+36other_tsrna.csv ./hc+36other_tsrna/median.csv
```


② 缺失值替换为NA
```r
setwd("D:/adult_dep/DE2")
data <- read.csv("./input_rename/hc+36other_tsrna.csv") 
data[data == 0] <- NA
write.csv(data, "./hc+36other_tsrna/tsrna_NA.csv", row.names = FALSE)  
```



③ Wilcoxon秩和检验（相关样本）  

秩和检验（Wilcoxon秩和检验）是一种用于比较两组相关或无关样本的非参数检验方法。它不要求数据符合正态分布，因此适用于那些无法满足正态性假设的情况。Wilcoxon秩和检验（相关样本）： 适用于配对设计，比如同一组个体在两个不同时间点的观测。在这种情况下，对每对相关样本进行秩和的比较。    
* 原理：在Wilcoxon秩和检验中，基本思想是将每个组的观测值按大小排序，然后为每个观测值分配一个秩次，最后通过比较秩和来判断两组样本是否来自同一分布。这个检验的零假设是两组样本来自相同的分布。  


```r
# 设置工作目录
setwd("D:/adult_dep/DE2")

# 读取数据
data <- read.csv("./hc+36other_tsrna/tsrna_NA.csv")

# 加载dplyr包
library(dplyr)

# 获取控制组和治疗组的列
control_columns <- grep("^control_", names(data), value = TRUE)
treatment_columns <- grep("^treatment_", names(data), value = TRUE)

# 执行 Wilcoxon 秩和检验
p_values_matrix <- matrix(NA, nrow = length(data$tsRNA), ncol = 3)
colnames(p_values_matrix) <- c("tsRNA", "p_value", "Q_value")

for (i in seq_along(data$tsRNA)) {
  tsRNA <- data$tsRNA[i]
  control_values <- unlist(data[data$tsRNA == tsRNA, control_columns])
  treatment_values <- unlist(data[data$tsRNA == tsRNA, treatment_columns])
  
  # 执行 Wilcoxon 秩和检验
  result <- try(wilcox.test(treatment_values, control_values))

  # 存储 p-value
  if (!inherits(result, "try-error")) {
    p_values_matrix[i, 1] <- tsRNA
    p_values_matrix[i, 2] <- result$p.value
  } else {
    cat("Error for tsRNA:", tsRNA, "\n")
  }
}

# 使用 Benjamini-Hochberg 方法进行多重比较校正
p_values_matrix[, 3] <- p.adjust(p_values_matrix[, 2], method = "BH")

# 输出结果
print(p_values_matrix)
# 将 p_values_matrix 保存到 CSV 文件
write.csv(p_values_matrix, file = "./hc+36other_tsrna/de_values_matrix.csv", row.names = FALSE)
```

④ 合并文件
```r
library(dplyr)

median_data <- read.csv("./hc+36other_tsrna/median.csv")
de_values_data <- read.csv("./hc+36other_tsrna/de_values_matrix.csv")
# 基于 miRNA 列合并数据
merged_data <- inner_join(median_data, de_values_data, by = "tsRNA")

# 查看合并后的数据
print(merged_data)
write.csv(merged_data, file = "./hc+36other_tsrna/result.csv", row.names = FALSE)

diff_gene <- subset(merged_data, Q_value < 0.05 & abs(log2fold_change) > 0.58)
write.csv(diff_gene, file="./hc+36other_tsrna/sig_q0.05_fc15.csv", quote = F)
```


⑤ 绘图
```r
# 加载包
library(ggplot2)
library(pheatmap)

merged_data$lgpadj <- -log10(merged_data$Q_value)

# 绘制火山图
volcano_plot <- ggplot(merged_data, aes(x = log2fold_change, y = lgpadj)) +
  geom_point(aes(color = ifelse(Q_value < 0.05 & log2fold_change > 1, "Up", 
                                ifelse(Q_value < 0.05 & log2fold_change < -1, "Down", "Not Sig"))), size = 2) +
  scale_color_manual(values = c("Up" = "#D32F2F", "Down" = "#1976D2", "Not Sig" = "grey70")) +
  labs(x = "Log2 Fold Change", y = "-Log10(padj)", title = "Volcano Plot") +
  theme_minimal(base_size = 14) +
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(),  # 去掉主网格线
        panel.grid.minor = element_blank(),  # 去掉次网格线
        axis.line = element_line(color = "black")) + # 显示横纵坐标轴
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", size = 0.5) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black", size = 0.5)

ggsave("./other_tsrna/volcano_plot.pdf", plot = volcano_plot, device = "pdf", dpi = 600, width = 8, height = 6)
```
```r
library(ggplot2)
library(pheatmap)

merged_data <- read.csv("./hc+36other_tsrna/result.csv")
merged_data$lgpadj <- -log10(merged_data$Q_value)

# Define categories for points
merged_data$category <- ifelse(merged_data$Q_value < 0.05 & merged_data$log2fold_change > 0.58, "Up", 
                                ifelse(merged_data$Q_value < 0.05 & merged_data$log2fold_change < -0.58, "Down", "Not Sig"))

# Calculate counts for each category
category_counts <- table(merged_data$category)

# Create a volcano plot
volcano_plot <- ggplot(merged_data, aes(x = log2fold_change, y = lgpadj)) +
  geom_point(aes(color = category), size = 2) +
  scale_color_manual(values = c("Up" = "#D32F2F", "Down" = "#1976D2", "Not Sig" = "grey70")) +
  labs(x = "Log2 Fold Change", y = "-Log10(padj)", title = "Volcano Plot") +
  theme_minimal(base_size = 14) +
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        axis.line = element_line(color = "black"),  # Show axis lines
        axis.ticks = element_line(color = "black"),  # Show axis ticks
        axis.ticks.length = unit(0.25, "cm")) +  # Set tick length
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", size = 0.5) +
  geom_vline(xintercept = c(-0.58, 0.58), linetype = "dashed", color = "black", size = 0.5) +
  
  # Add category count labels
  geom_text(aes(x = 1, y = 30, label = paste("Up:", category_counts["Up"])), color = "#D32F2F", size = 5, hjust = 0) +
  geom_text(aes(x = 1, y = 28, label = paste("Down:", category_counts["Down"])), color = "#1976D2", size = 5, hjust = 0) +
  geom_text(aes(x = 1, y = 26, label = paste("Not Sig:", category_counts["Not Sig"])), color = "grey70", size = 5, hjust = 0)

# Save the plot
ggsave("./hc+36other_tsrna/hc+36other_tsrna_volcano_plot.pdf", plot = volcano_plot, device = "pdf", dpi = 600, width = 8, height = 6)
```




* 图中标注

```r

# 图中标注名称
annotations <- data.frame(
  tsRNA = c("TTGGTCGTGGTTGTAGTCCGTGCGAGAATACCA", "GTTTCCGTAGTGTAGTGGTTATCACGTT", "GCGGGAGACCGGGGTTCGATTCCCCGACGGGGAGC", "CATTGGTCGTGGTTGTAGTCCGTGCGAGAATACCA", "CCCTGGTGGTCTAGTGGTTAGGATTCGGCGCT"),
  new_label = c("tsRNA-#98", "tsRNA-#57", "tsRNA-#51", "tsRNA-#23", "tsRNA-#25")
)

# 提取需要标注的 miRNA 数据
annotations <- merged_data[merged_data$tsRNA %in% annotations$tsRNA, ]

# Add annotations for the top 5 miRNAs with the highest p-value (lowest significance)
top_p_values <- merged_data[order(merged_data$lgpadj, decreasing = TRUE), ]
top_p_annotations <- head(top_p_values, 3)

# Combine the annotations: predefined list and the top 5 p-value miRNAs
all_annotations <- rbind(annotations, top_p_annotations)

# Create volcano plot
volcano_plot <- ggplot(merged_data, aes(x = log2fold_change, y = lgpadj)) +
  geom_point(aes(color = ifelse(Q_value < 0.05 & log2fold_change > 0.58, "Up", 
                                ifelse(Q_value < 0.05 & log2fold_change < -0.58, "Down", "Not Sig"))), size = 2) +
  scale_color_manual(values = c("Up" = "#D32F2F", "Down" = "#1976D2", "Not Sig" = "grey70")) +
  labs(x = "Log2 Fold Change", y = "-Log10(padj)", title = "Volcano Plot") +
  theme_minimal(base_size = 14) +
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", size = 0.5) +
  geom_vline(xintercept = c(-0.58, 0.58), linetype = "dashed", color = "black", size = 0.5) +
  geom_point(data = annotations, aes(color = "Highlighted"), size = 3, shape = 1, stroke = 1.5, fill = NA) +  # Highlight miRNAs from predefined list
  geom_text(data = all_annotations, aes(label = tsRNA), vjust = -1, hjust = 0.5, size = 3, color = "black")  # Add text labels for both predefined and top p-value miRNAs

# Display the volcano plot
print(volcano_plot)


ggsave("./hc+36other_tsrna/label_hc+36other_tsrna_volcano_plot3.pdf", plot = volcano_plot, device = "pdf", dpi = 600, width = 16, height = 12)

ggsave("./hc+36other_tsrna/label_hc+36other_tsrna_volcano_plot2.pdf", plot = volcano_plot, device = "pdf", dpi = 600, width = 8, height = 6)
```

```r
# 绘制热图
library(ggplot2)
library(tidyr)
library(dplyr)

# 读取数据
merged_data <- read.csv("./other_tsrna/result.csv")

# 指定需要标注的 tsRNA 和对应的新名称
tsRNAs_to_annotate <- c("ATCCTGCCGACTACGCCA", "TCCCTGGTGGTCTAGTGGTTAGGATTCGGCGCT", "GCGGGAGACCGGGGTTCGATTCCCCGACGGGGAGC")

new_names <- c("tsRNA-#17", "tsRNA-#83", "tsRNA-#51")  # 新名称

# 创建一个命名映射
name_mapping <- setNames(new_names, tsRNAs_to_annotate)

# 提取需要的 tsRNA 数据并设置顺序
heatmap_data <- merged_data %>%
  filter(tsRNA %in% tsRNAs_to_annotate) %>%
  select(tsRNA, log2control_median, log2treatment_median) %>%
  pivot_longer(cols = starts_with("log2"), names_to = "Condition", values_to = "Expression") %>%
  mutate(tsRNA = factor(name_mapping[tsRNA], levels = new_names))  # 使用新名称

# 定义颜色从白色到深蓝色
start_color <- rgb(255/255, 255/255, 255/255)  # 白色
end_color <- rgb(0/255, 51/255, 102/255)        # 深蓝色

# 生成渐变
num_colors <- 10
gradient_colors <- colorRampPalette(c(start_color, end_color))(num_colors)

# 绘制热图
heatmap_plot <- ggplot(heatmap_data, aes(x = tsRNA, y = Condition, fill = Expression)) +
  geom_tile(color = "white") +
  scale_fill_gradientn(colors = gradient_colors, 
                       name = "Expression Level") +
  labs(x = "tsRNA", y = "Condition", title = "Heatmap of tsRNA Expression") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 保存热图
ggsave("./other_tsrna/5tsRNA_heatmap.pdf", plot = heatmap_plot, device = "pdf", dpi = 600, width = 8, height = 6)
```


































## 1. 36other_mirna

2. 标注control_vs_treatment，手动标注
```r
setwd("D:/adult_dep/DE2")
df <- read.csv("./input/rpm_50%_log_mirna.csv", header = TRUE, row.names = 1)
df <- df[, df[1, ] %in% c(1, 3)]
colnames(df) <- ifelse(df[1, ] == 3, 
                       paste0("control_", colnames(df)), 
                       paste0("treatment_", colnames(df)))
control_columns <- colnames(df)[grepl("^control_", colnames(df))]
treatment_columns <- colnames(df)[grepl("^treatment_", colnames(df))]        

# 删除第一行
df <- df[-1, ]
# 查看修改后的数据框
head(df)
write.csv(df, file = "./input_rename/36other_mirna.csv")
```

3. 找到每个mirna的treatment_vs_control中位数，计算FC
```bash
cd /mnt/d/adult_dep/DE2
# 记得加“mirna”
python median_mirna.py ./input_rename/36other_mirna.csv ./36other_mirna/median.csv
```


② 缺失值替换为NA
```r
setwd("D:/adult_dep/DE2")
data <- read.csv("./input_rename/36other_mirna.csv") 
data[data == 0] <- NA
write.csv(data, "./36other_mirna/mirna_NA.csv", row.names = FALSE)  
```



③ Wilcoxon秩和检验（相关样本）  

秩和检验（Wilcoxon秩和检验）是一种用于比较两组相关或无关样本的非参数检验方法。它不要求数据符合正态分布，因此适用于那些无法满足正态性假设的情况。Wilcoxon秩和检验（相关样本）： 适用于配对设计，比如同一组个体在两个不同时间点的观测。在这种情况下，对每对相关样本进行秩和的比较。    
* 原理：在Wilcoxon秩和检验中，基本思想是将每个组的观测值按大小排序，然后为每个观测值分配一个秩次，最后通过比较秩和来判断两组样本是否来自同一分布。这个检验的零假设是两组样本来自相同的分布。  


```r
# 设置工作目录
setwd("D:/adult_dep/DE2")

# 读取数据
data <- read.csv("./36other_mirna/mirna_NA.csv")

# 加载dplyr包
library(dplyr)

# 获取控制组和治疗组的列
control_columns <- grep("^control_", names(data), value = TRUE)
treatment_columns <- grep("^treatment_", names(data), value = TRUE)


# 执行 Wilcoxon 秩和检验
p_values_matrix <- matrix(NA, nrow = length(data$miRNA), ncol = 3)
colnames(p_values_matrix) <- c("miRNA", "p_value", "Q_value")

for (i in seq_along(data$miRNA)) {
  miRNA <- data$miRNA[i]
  control_values <- unlist(data[data$miRNA == miRNA, control_columns])
  treatment_values <- unlist(data[data$miRNA == miRNA, treatment_columns])
  
  # 执行 Wilcoxon 秩和检验
  result <- try(wilcox.test(treatment_values, control_values))

  # 存储 p-value
  if (!inherits(result, "try-error")) {
    p_values_matrix[i, 1] <- miRNA
    p_values_matrix[i, 2] <- result$p.value
  } else {
    cat("Error for miRNA:", miRNA, "\n")
  }
}

# 使用 Benjamini-Hochberg 方法进行多重比较校正
p_values_matrix[, 3] <- p.adjust(p_values_matrix[, 2], method = "BH")

# 输出结果
print(p_values_matrix)
# 将 p_values_matrix 保存到 CSV 文件
write.csv(p_values_matrix, file = "./36other_mirna/de_values_matrix.csv", row.names = FALSE)
```

④ 合并文件
```r
library(dplyr)

median_data <- read.csv("./36other_mirna/median.csv")
de_values_data <- read.csv("./36other_mirna/de_values_matrix.csv")
# 基于 miRNA 列合并数据
merged_data <- inner_join(median_data, de_values_data, by = "miRNA")

# 查看合并后的数据
print(merged_data)
write.csv(merged_data, file = "./36other_mirna/result.csv", row.names = FALSE)

diff_gene <- subset(merged_data, Q_value < 0.05 & abs(log2fold_change) > 0.58)
write.csv(diff_gene, file="./36other_mirna/sig_q0.05_fc15.csv", quote = F)
```









⑤ 绘图
```r
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}
if (!requireNamespace("pheatmap", quietly = TRUE)) {
  install.packages("pheatmap")
}

# 加载包
library(ggplot2)
library(pheatmap)

merged_data$lgpadj <- -log10(merged_data$Q_value)

# 绘制火山图
volcano_plot <- ggplot(merged_data, aes(x = log2fold_change, y = lgpadj)) +
  geom_point(aes(color = ifelse(Q_value < 0.05 & log2fold_change > 1, "Up", 
                                ifelse(Q_value < 0.05 & log2fold_change < -1, "Down", "Not Sig"))), size = 2) +
  scale_color_manual(values = c("Up" = "#D32F2F", "Down" = "#1976D2", "Not Sig" = "grey70")) +
  labs(x = "Log2 Fold Change", y = "-Log10(padj)", title = "Volcano Plot") +
  theme_minimal(base_size = 14) +
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(),  # 去掉主网格线
        panel.grid.minor = element_blank(),  # 去掉次网格线
        axis.line = element_line(color = "black")) + # 显示横纵坐标轴
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", size = 0.5) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black", size = 0.5)

ggsave("./bp_mirna/volcano_plot.pdf", plot = volcano_plot, device = "pdf", dpi = 600, width = 8, height = 6)
```

```r
# 加载包
library(ggplot2)
library(pheatmap)

merged_data <- read.csv("./36other_mirna/result.csv")
merged_data$lgpadj <- -log10(merged_data$Q_value)

# Define categories for points
merged_data$category <- ifelse(merged_data$Q_value < 0.05 & merged_data$log2fold_change > 0.58, "Up", 
                                ifelse(merged_data$Q_value < 0.05 & merged_data$log2fold_change < -0.58, "Down", "Not Sig"))

# Calculate counts for each category
category_counts <- table(merged_data$category)

# Create a volcano plot
volcano_plot <- ggplot(merged_data, aes(x = log2fold_change, y = lgpadj)) +
  geom_point(aes(color = category), size = 2) +
  scale_color_manual(values = c("Up" = "#D32F2F", "Down" = "#1976D2", "Not Sig" = "grey70")) +
  labs(x = "Log2 Fold Change", y = "-Log10(padj)", title = "Volcano Plot") +
  theme_minimal(base_size = 14) +
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        axis.line = element_line(color = "black"),  # Show axis lines
        axis.ticks = element_line(color = "black"),  # Show axis ticks
        axis.ticks.length = unit(0.25, "cm")) +  # Set tick length
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", size = 0.5) +
  geom_vline(xintercept = c(-0.58, 0.58), linetype = "dashed", color = "black", size = 0.5) +
  
  # Add category count labels
  geom_text(aes(x = 1, y = 30, label = paste("Up:", category_counts["Up"])), color = "#D32F2F", size = 5, hjust = 0) +
  geom_text(aes(x = 1, y = 28, label = paste("Down:", category_counts["Down"])), color = "#1976D2", size = 5, hjust = 0) +
  geom_text(aes(x = 1, y = 26, label = paste("Not Sig:", category_counts["Not Sig"])), color = "grey70", size = 5, hjust = 0)

# Save the plot
ggsave("./36other_mirna/36other_mirna_volcano_plot.pdf", plot = volcano_plot, device = "pdf", dpi = 600, width = 8, height = 6)
```

```r
# 图中标注名称
miRNAs_to_annotate <- c("hsa-mir-28_hsa-miR-28-3p", "hsa-mir-410_hsa-miR-410-3p", 
                         "hsa-mir-4433b_hsa-miR-4433b-5p")

# 提取需要标注的 miRNA 数据
annotations <- merged_data[merged_data$miRNA %in% miRNAs_to_annotate, ]

# 创建火山图
volcano_plot <- ggplot(merged_data, aes(x = log2fold_change, y = lgpadj)) +
  geom_point(aes(color = ifelse(Q_value < 0.05 & log2fold_change > 0.58, "Up", 
                                ifelse(Q_value < 0.05 & log2fold_change < -0.58, "Down", "Not Sig"))), size = 2) +
  scale_color_manual(values = c("Up" = "#D32F2F", "Down" = "#1976D2", "Not Sig" = "grey70")) +
  labs(x = "Log2 Fold Change", y = "-Log10(padj)", title = "Volcano Plot") +
  theme_minimal(base_size = 14) +
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", size = 0.5) +
  geom_vline(xintercept = c(-0.58, 0.58), linetype = "dashed", color = "black", size = 0.5) +
  geom_point(data = annotations, aes(color = "Highlighted"), size = 3, shape = 1, stroke = 1.5, fill = NA) +  # 添加圈圈
  geom_text(data = annotations, aes(label = miRNA), vjust = -1, hjust = 0.5, size = 3, color = "black")


ggsave("./bp_mirna/label_bp_mirna_volcano_plot.pdf", plot = volcano_plot, device = "pdf", dpi = 600, width = 8, height = 6)

```

```r
# 绘制热图
library(ggplot2)
library(tidyr)
library(dplyr)

# 读取数据
merged_data <- read.csv("./bp_mirna/result.csv")

# 指定要绘制热图的 miRNA
miRNAs_to_annotate <- c("hsa-mir-28_hsa-miR-28-3p", "hsa-mir-410_hsa-miR-410-3p", 
                         "hsa-mir-4433b_hsa-miR-4433b-5p")

# 提取需要的 miRNA 数据并设置顺序
heatmap_data <- merged_data %>%
  filter(miRNA %in% miRNAs_to_annotate) %>%
  select(miRNA, log2control_median, log2treatment_median) %>%
  pivot_longer(cols = starts_with("log2"), names_to = "Condition", values_to = "Expression") %>%
  mutate(miRNA = factor(miRNA, levels = miRNAs_to_annotate))  # 设置因子顺序

# 定义颜色从白色到深蓝色
start_color <- rgb(255/255, 255/255, 255/255)  # 白色
end_color <- rgb(0/255, 51/255, 102/255)        # 深蓝色

# 生成渐变
num_colors <- 10
gradient_colors <- colorRampPalette(c(start_color, end_color))(num_colors)

# 绘制热图
heatmap_plot <- ggplot(heatmap_data, aes(x = miRNA, y = Condition, fill = Expression)) +
  geom_tile(color = "white") +
  scale_fill_gradientn(colors = gradient_colors, 
                       name = "Expression Level") +
  labs(x = "miRNA", y = "Condition", title = "Heatmap of miRNA Expression") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 保存热图
ggsave("./bp_mirna/5mirna_heatmap.pdf", plot = heatmap_plot, device = "pdf", dpi = 600, width = 8, height = 6)
```





### 2. 36other_pirna

2. 标注control_vs_treatment，手动标注
```r
setwd("D:/adult_dep/DE2")
df <- read.csv("./input/rpm_50%_log_pirna.csv", header = TRUE, row.names = 1)
df <- df[, df[1, ] %in% c(1, 3)]
colnames(df) <- ifelse(df[1, ] == 3, 
                       paste0("control_", colnames(df)), 
                       paste0("treatment_", colnames(df)))
control_columns <- colnames(df)[grepl("^control_", colnames(df))]        

# 删除第一行
df <- df[-1, ]
# 查看修改后的数据框
head(df)
write.csv(df, file = "./input_rename/36other_pirna.csv")
```

3. 找到每个mirna的treatment_vs_control中位数，计算FC
```bash
cd /mnt/d/adult_dep/DE2
# 记得加“pirna”
python median_pirna.py ./input_rename/36other_pirna.csv ./36other_pirna/median.csv
```


② 缺失值替换为NA
```r
setwd("D:/adult_dep/DE2")
data <- read.csv("./input_rename/36other_pirna.csv") 
data[data == 0] <- NA
write.csv(data, "./36other_pirna/pirna_NA.csv", row.names = FALSE)  
```



③ Wilcoxon秩和检验（相关样本）  

秩和检验（Wilcoxon秩和检验）是一种用于比较两组相关或无关样本的非参数检验方法。它不要求数据符合正态分布，因此适用于那些无法满足正态性假设的情况。Wilcoxon秩和检验（相关样本）： 适用于配对设计，比如同一组个体在两个不同时间点的观测。在这种情况下，对每对相关样本进行秩和的比较。    
* 原理：在Wilcoxon秩和检验中，基本思想是将每个组的观测值按大小排序，然后为每个观测值分配一个秩次，最后通过比较秩和来判断两组样本是否来自同一分布。这个检验的零假设是两组样本来自相同的分布。  


```r
# 设置工作目录
setwd("D:/adult_dep/DE2")

# 读取数据
data <- read.csv("./36other_pirna/pirna_NA.csv")

# 加载dplyr包
library(dplyr)

# 获取控制组和治疗组的列
control_columns <- grep("^control_", names(data), value = TRUE)
treatment_columns <- grep("^treatment_", names(data), value = TRUE)
# 835

# 执行 Wilcoxon 秩和检验
p_values_matrix <- matrix(NA, nrow = length(data$piRNA), ncol = 3)
colnames(p_values_matrix) <- c("piRNA", "p_value", "Q_value")

for (i in seq_along(data$piRNA)) {
  piRNA <- data$piRNA[i]
  control_values <- unlist(data[data$piRNA == piRNA, control_columns])
  treatment_values <- unlist(data[data$piRNA == piRNA, treatment_columns])
  
  # 执行 Wilcoxon 秩和检验
  result <- try(wilcox.test(treatment_values, control_values))

  # 存储 p-value
  if (!inherits(result, "try-error")) {
    p_values_matrix[i, 1] <- piRNA
    p_values_matrix[i, 2] <- result$p.value
  } else {
    cat("Error for piRNA:", piRNA, "\n")
  }
}

# 使用 Benjamini-Hochberg 方法进行多重比较校正
p_values_matrix[, 3] <- p.adjust(p_values_matrix[, 2], method = "BH")

# 输出结果
print(p_values_matrix)
# 将 p_values_matrix 保存到 CSV 文件
write.csv(p_values_matrix, file = "./36other_pirna/de_values_matrix.csv", row.names = FALSE)
```

④ 合并文件
```r
library(dplyr)

median_data <- read.csv("./36other_pirna/median.csv")
de_values_data <- read.csv("./36other_pirna/de_values_matrix.csv")
# 基于 miRNA 列合并数据
merged_data <- inner_join(median_data, de_values_data, by = "piRNA")

# 查看合并后的数据
print(merged_data)
write.csv(merged_data, file = "./36other_pirna/result.csv", row.names = FALSE)

diff_gene <- subset(merged_data, Q_value < 0.05 & abs(log2fold_change) > 0.58)
write.csv(diff_gene, file="./36other_pirna/sig_q0.05_fc15.csv", quote = F)
```




⑤ 绘图
```r
# 加载包
library(ggplot2)
library(pheatmap)

merged_data$lgpadj <- -log10(merged_data$Q_value)

# 绘制火山图
volcano_plot <- ggplot(merged_data, aes(x = log2fold_change, y = lgpadj)) +
  geom_point(aes(color = ifelse(Q_value < 0.05 & log2fold_change > 1, "Up", 
                                ifelse(Q_value < 0.05 & log2fold_change < -1, "Down", "Not Sig"))), size = 2) +
  scale_color_manual(values = c("Up" = "#D32F2F", "Down" = "#1976D2", "Not Sig" = "grey70")) +
  labs(x = "Log2 Fold Change", y = "-Log10(padj)", title = "Volcano Plot") +
  theme_minimal(base_size = 14) +
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(),  # 去掉主网格线
        panel.grid.minor = element_blank(),  # 去掉次网格线
        axis.line = element_line(color = "black")) + # 显示横纵坐标轴
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", size = 0.5) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black", size = 0.5)

ggsave("./bp_pirna/volcano_plot.pdf", plot = volcano_plot, device = "pdf", dpi = 600, width = 8, height = 6)
```




```r
library(ggplot2)
library(pheatmap)

merged_data <- read.csv("./36other_pirna/result.csv")
merged_data$lgpadj <- -log10(merged_data$Q_value)

# Define categories for points
merged_data$category <- ifelse(merged_data$Q_value < 0.05 & merged_data$log2fold_change > 0.58, "Up", 
                                ifelse(merged_data$Q_value < 0.05 & merged_data$log2fold_change < -0.58, "Down", "Not Sig"))

# Calculate counts for each category
category_counts <- table(merged_data$category)

# Create a volcano plot
volcano_plot <- ggplot(merged_data, aes(x = log2fold_change, y = lgpadj)) +
  geom_point(aes(color = category), size = 2) +
  scale_color_manual(values = c("Up" = "#D32F2F", "Down" = "#1976D2", "Not Sig" = "grey70")) +
  labs(x = "Log2 Fold Change", y = "-Log10(padj)", title = "Volcano Plot") +
  theme_minimal(base_size = 14) +
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        axis.line = element_line(color = "black"),  # Show axis lines
        axis.ticks = element_line(color = "black"),  # Show axis ticks
        axis.ticks.length = unit(0.25, "cm")) +  # Set tick length
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", size = 0.5) +
  geom_vline(xintercept = c(-0.58, 0.58), linetype = "dashed", color = "black", size = 0.5) +
  
  # Add category count labels
  geom_text(aes(x = 1, y = 30, label = paste("Up:", category_counts["Up"])), color = "#D32F2F", size = 5, hjust = 0) +
  geom_text(aes(x = 1, y = 28, label = paste("Down:", category_counts["Down"])), color = "#1976D2", size = 5, hjust = 0) +
  geom_text(aes(x = 1, y = 26, label = paste("Not Sig:", category_counts["Not Sig"])), color = "grey70", size = 5, hjust = 0)

# Save the plot
ggsave("./36other_pirna/36other_pirna_volcano_plot.pdf", plot = volcano_plot, device = "pdf", dpi = 600, width = 8, height = 6)
```
* 图中标注

```r
# 图中标注名称
miRNAs_to_annotate <- c("piR-hsa-135526", "piR-hsa-147698","piR-hsa-151782","piR-hsa-184607","piR-hsa-789156")

# 提取需要标注的 miRNA 数据
annotations <- merged_data[merged_data$piRNA %in% miRNAs_to_annotate, ]

# 创建火山图
volcano_plot <- ggplot(merged_data, aes(x = log2fold_change, y = lgpadj)) +
  geom_point(aes(color = ifelse(Q_value < 0.05 & log2fold_change > 0.58, "Up", 
                                ifelse(Q_value < 0.05 & log2fold_change < -0.58, "Down", "Not Sig"))), size = 2) +
  scale_color_manual(values = c("Up" = "#D32F2F", "Down" = "#1976D2", "Not Sig" = "grey70")) +
  labs(x = "Log2 Fold Change", y = "-Log10(padj)", title = "Volcano Plot") +
  theme_minimal(base_size = 14) +
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", size = 0.5) +
  geom_vline(xintercept = c(-0.58, 0.58), linetype = "dashed", color = "black", size = 0.5) +
  geom_point(data = annotations, aes(color = "Highlighted"), size = 3, shape = 1, stroke = 1.5, fill = NA) +  # 添加圈圈
  geom_text(data = annotations, aes(label = piRNA), vjust = -1, hjust = 0.5, size = 3, color = "black")


ggsave("./bp_pirna/label_bp_pirna_volcano_plot.pdf", plot = volcano_plot, device = "pdf", dpi = 600, width = 8, height = 6)

```

```r
# 绘制热图
library(ggplot2)
library(tidyr)
library(dplyr)

# 读取数据
merged_data <- read.csv("./bp_pirna/result.csv")

# 指定要绘制热图的 miRNA
piRNAs_to_annotate <- c("piR-hsa-135526", "piR-hsa-147698","piR-hsa-151782","piR-hsa-184607","piR-hsa-789156")

# 提取需要的 piRNA 数据并设置顺序
heatmap_data <- merged_data %>%
  filter(piRNA %in% piRNAs_to_annotate) %>%
  select(piRNA, log2control_median, log2treatment_median) %>%
  pivot_longer(cols = starts_with("log2"), names_to = "Condition", values_to = "Expression") %>%
  mutate(piRNA = factor(piRNA, levels = piRNAs_to_annotate))  # 设置因子顺序

# 定义颜色从白色到深蓝色
start_color <- rgb(255/255, 255/255, 255/255)  # 白色
end_color <- rgb(0/255, 51/255, 102/255)        # 深蓝色

# 生成渐变
num_colors <- 10
gradient_colors <- colorRampPalette(c(start_color, end_color))(num_colors)

# 绘制热图
heatmap_plot <- ggplot(heatmap_data, aes(x = piRNA, y = Condition, fill = Expression)) +
  geom_tile(color = "white") +
  scale_fill_gradientn(colors = gradient_colors, 
                       name = "Expression Level") +
  labs(x = "piRNA", y = "Condition", title = "Heatmap of piRNA Expression") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 保存热图
ggsave("./bp_pirna/5piRNA_heatmap.pdf", plot = heatmap_plot, device = "pdf", dpi = 600, width = 8, height = 6)
```









## 3. 36other_rsrna

2. 标注control_vs_treatment，手动标注
```r
setwd("D:/adult_dep/DE2")
df <- read.csv("./input/rpm_50%_log_rsrna.csv", header = TRUE, row.names = 1)
df <- df[, df[1, ] %in% c(1, 3)]
colnames(df) <- ifelse(df[1, ] == 3, 
                       paste0("control_", colnames(df)), 
                       paste0("treatment_", colnames(df)))
control_columns <- colnames(df)[grepl("^control_", colnames(df))]        

# 删除第一行
df <- df[-1, ]
# 查看修改后的数据框
head(df)
write.csv(df, file = "./input_rename/36other_rsrna.csv")
```

3. 找到每个mirna的treatment_vs_control中位数，计算FC
```bash
cd /mnt/d/adult_dep/DE2
# 记得加“rsrna”
python median_rsrna.py ./input_rename/36other_rsrna.csv ./36other_rsrna/median.csv
```


② 缺失值替换为NA
```r
setwd("D:/adult_dep/DE2")
data <- read.csv("./input_rename/36other_rsrna.csv") 
data[data == 0] <- NA
write.csv(data, "./36other_rsrna/rsrna_NA.csv", row.names = FALSE)  
```



③ Wilcoxon秩和检验（相关样本）  

秩和检验（Wilcoxon秩和检验）是一种用于比较两组相关或无关样本的非参数检验方法。它不要求数据符合正态分布，因此适用于那些无法满足正态性假设的情况。Wilcoxon秩和检验（相关样本）： 适用于配对设计，比如同一组个体在两个不同时间点的观测。在这种情况下，对每对相关样本进行秩和的比较。    
* 原理：在Wilcoxon秩和检验中，基本思想是将每个组的观测值按大小排序，然后为每个观测值分配一个秩次，最后通过比较秩和来判断两组样本是否来自同一分布。这个检验的零假设是两组样本来自相同的分布。  


```r
# 设置工作目录
setwd("D:/adult_dep/DE2")

# 读取数据
data <- read.csv("./36other_rsrna/rsrna_NA.csv")

# 加载dplyr包
library(dplyr)

# 获取控制组和治疗组的列
control_columns <- grep("^control_", names(data), value = TRUE)
treatment_columns <- grep("^treatment_", names(data), value = TRUE)

# 执行 Wilcoxon 秩和检验
p_values_matrix <- matrix(NA, nrow = length(data$rsRNA), ncol = 3)
colnames(p_values_matrix) <- c("rsRNA", "p_value", "Q_value")

for (i in seq_along(data$rsRNA)) {
  rsRNA <- data$rsRNA[i]
  control_values <- unlist(data[data$rsRNA == rsRNA, control_columns])
  treatment_values <- unlist(data[data$rsRNA == rsRNA, treatment_columns])
  
  # 执行 Wilcoxon 秩和检验
  result <- try(wilcox.test(treatment_values, control_values))

  # 存储 p-value
  if (!inherits(result, "try-error")) {
    p_values_matrix[i, 1] <- rsRNA
    p_values_matrix[i, 2] <- result$p.value
  } else {
    cat("Error for rsRNA:", rsRNA, "\n")
  }
}

# 使用 Benjamini-Hochberg 方法进行多重比较校正
p_values_matrix[, 3] <- p.adjust(p_values_matrix[, 2], method = "BH")

# 输出结果
print(p_values_matrix)
write.csv(p_values_matrix, file = "./36other_rsrna/de_values_matrix.csv", row.names = FALSE)
```

④ 合并文件
```r
library(dplyr)

median_data <- read.csv("./36other_rsrna/median.csv")
de_values_data <- read.csv("./36other_rsrna/de_values_matrix.csv")
# 基于 miRNA 列合并数据
merged_data <- inner_join(median_data, de_values_data, by = "rsRNA")

# 查看合并后的数据
print(merged_data)
write.csv(merged_data, file = "./36other_rsrna/result.csv", row.names = FALSE)

diff_gene <- subset(merged_data, Q_value < 0.05 & abs(log2fold_change) > 0.58)
write.csv(diff_gene, file="./36other_rsrna/sig_q0.05_fc15.csv", quote = F)
```


⑤ 绘图
```r
# 加载包
library(ggplot2)
library(pheatmap)

merged_data$lgpadj <- -log10(merged_data$Q_value)

# 绘制火山图
volcano_plot <- ggplot(merged_data, aes(x = log2fold_change, y = lgpadj)) +
  geom_point(aes(color = ifelse(Q_value < 0.05 & log2fold_change > 1, "Up", 
                                ifelse(Q_value < 0.05 & log2fold_change < -1, "Down", "Not Sig"))), size = 2) +
  scale_color_manual(values = c("Up" = "#D32F2F", "Down" = "#1976D2", "Not Sig" = "grey70")) +
  labs(x = "Log2 Fold Change", y = "-Log10(padj)", title = "Volcano Plot") +
  theme_minimal(base_size = 14) +
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(),  # 去掉主网格线
        panel.grid.minor = element_blank(),  # 去掉次网格线
        axis.line = element_line(color = "black")) + # 显示横纵坐标轴
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", size = 0.5) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black", size = 0.5)

ggsave("./bp_rsrna/volcano_plot.pdf", plot = volcano_plot, device = "pdf", dpi = 600, width = 8, height = 6)
```

```r
library(ggplot2)
library(pheatmap)

merged_data <- read.csv("./36other_rsrna/result.csv")
merged_data$lgpadj <- -log10(merged_data$Q_value)

# Define categories for points
merged_data$category <- ifelse(merged_data$Q_value < 0.05 & merged_data$log2fold_change > 0.58, "Up", 
                                ifelse(merged_data$Q_value < 0.05 & merged_data$log2fold_change < -0.58, "Down", "Not Sig"))

# Calculate counts for each category
category_counts <- table(merged_data$category)

# Create a volcano plot
volcano_plot <- ggplot(merged_data, aes(x = log2fold_change, y = lgpadj)) +
  geom_point(aes(color = category), size = 2) +
  scale_color_manual(values = c("Up" = "#D32F2F", "Down" = "#1976D2", "Not Sig" = "grey70")) +
  labs(x = "Log2 Fold Change", y = "-Log10(padj)", title = "Volcano Plot") +
  theme_minimal(base_size = 14) +
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        axis.line = element_line(color = "black"),  # Show axis lines
        axis.ticks = element_line(color = "black"),  # Show axis ticks
        axis.ticks.length = unit(0.25, "cm")) +  # Set tick length
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", size = 0.5) +
  geom_vline(xintercept = c(-0.58, 0.58), linetype = "dashed", color = "black", size = 0.5) +
  
  # Add category count labels
  geom_text(aes(x = 1, y = 30, label = paste("Up:", category_counts["Up"])), color = "#D32F2F", size = 5, hjust = 0) +
  geom_text(aes(x = 1, y = 28, label = paste("Down:", category_counts["Down"])), color = "#1976D2", size = 5, hjust = 0) +
  geom_text(aes(x = 1, y = 26, label = paste("Not Sig:", category_counts["Not Sig"])), color = "grey70", size = 5, hjust = 0)

# Save the plot
ggsave("./36other_rsrna/36other_rsrna_volcano_plot.pdf", plot = volcano_plot, device = "pdf", dpi = 600, width = 8, height = 6)
```




* 图中标注

```r
# 图中标注名称
annotations <- data.frame(
  rsRNA = c("ACAATGTAGGTAAGGGAAGTCG", "ACGGACCAAGGAGTCTAACACGTGCGCG", "CGACTCTTAGCGGTGGATCACTCGGCTCGTGCG","GGAACAATGTAGGTAAGGGA","GGAGTCTAACACGTGCGCG","TGGGAGACCGCCTGGGAATACCGGGTGCTGTAGGCTT", "GACTCTTAGCGGTGGATCACTCGGCTCGTG"),
  new_label = c("rsRNA-#18", "rsRNA-#29", "rsRNA-#95","rsRNA-#176","rsRNA-#182","rsRNA-#261", "rsRNA-#155")
)

# 提取需要标注的 miRNA 数据
annotations_data <- merged_data[merged_data$rsRNA %in% annotations$rsRNA, ]

# 创建火山图
volcano_plot <- ggplot(merged_data, aes(x = log2fold_change, y = lgpadj)) +
  geom_point(aes(color = ifelse(Q_value < 0.05 & log2fold_change > 0.58, "Up", 
                                ifelse(Q_value < 0.05 & log2fold_change < -0.58, "Down", "Not Sig"))), size = 2) +
  scale_color_manual(values = c("Up" = "#D32F2F", "Down" = "#1976D2", "Not Sig" = "grey70")) +
  labs(x = "Log2 Fold Change", y = "-Log10(padj)", title = "Volcano Plot") +
  theme_minimal(base_size = 14) +
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", size = 0.5) +
  geom_vline(xintercept = c(-0.58, 0.58), linetype = "dashed", color = "black", size = 0.5) +
  geom_point(data = annotations_data, aes(color = "Highlighted"), size = 3, shape = 1, stroke = 1.5, fill = NA) +  # 添加圈圈
  geom_text(data = merge(annotations_data, annotations, by = "rsRNA"), aes(label = new_label), 
            vjust = -1, hjust = 0.5, size = 3, color = "black")  # 使用新的标签


ggsave("./bp_rsrna/label_bp_rsrna_volcano_plot.pdf", plot = volcano_plot, device = "pdf", dpi = 600, width = 8, height = 6)

```

```r
# 绘制热图
library(ggplot2)
library(tidyr)
library(dplyr)

# 读取数据
merged_data <- read.csv("./bp_rsrna/result.csv")

# 指定需要标注的 tsRNA 和对应的新名称
rsRNAs_to_annotate <- c("ACAATGTAGGTAAGGGAAGTCG", "ACGGACCAAGGAGTCTAACACGTGCGCG", "CGACTCTTAGCGGTGGATCACTCGGCTCGTGCG","GGAACAATGTAGGTAAGGGA","GGAGTCTAACACGTGCGCG","TGGGAGACCGCCTGGGAATACCGGGTGCTGTAGGCTT", "GACTCTTAGCGGTGGATCACTCGGCTCGTG")

new_names <- c("rsRNA-#18", "rsRNA-#29", "rsRNA-#95","rsRNA-#176","rsRNA-#182","rsRNA-#261", "rsRNA-#155")

# 创建一个命名映射
name_mapping <- setNames(new_names, rsRNAs_to_annotate)

# 提取需要的 rsRNA 数据并设置顺序
heatmap_data <- merged_data %>%
  filter(rsRNA %in% rsRNAs_to_annotate) %>%
  select(rsRNA, log2control_median, log2treatment_median) %>%
  pivot_longer(cols = starts_with("log2"), names_to = "Condition", values_to = "Expression") %>%
  mutate(rsRNA = factor(name_mapping[rsRNA], levels = new_names))  # 使用新名称

# 定义颜色从白色到深蓝色
start_color <- rgb(255/255, 255/255, 255/255)  # 白色
end_color <- rgb(0/255, 51/255, 102/255)        # 深蓝色

# 生成渐变
num_colors <- 10
gradient_colors <- colorRampPalette(c(start_color, end_color))(num_colors)

# 绘制热图
heatmap_plot <- ggplot(heatmap_data, aes(x = rsRNA, y = Condition, fill = Expression)) +
  geom_tile(color = "white") +
  scale_fill_gradientn(colors = gradient_colors, 
                       name = "Expression Level") +
  labs(x = "rsRNA", y = "Condition", title = "Heatmap of rsRNA Expression") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 保存热图
ggsave("./bp_rsrna/5bp_heatmap.pdf", plot = heatmap_plot, device = "pdf", dpi = 600, width = 8, height = 6)
```





## 4. 36other_tsrna

2. 标注control_vs_treatment，手动标注
```r
setwd("D:/adult_dep/DE2")
df <- read.csv("./input/rpm_50%_log_tsrna.csv", header = TRUE, row.names = 1)
df <- df[, df[1, ] %in% c(1, 3)]
colnames(df) <- ifelse(df[1, ] == 3, 
                       paste0("control_", colnames(df)), 
                       paste0("treatment_", colnames(df)))
control_columns <- colnames(df)[grepl("^control_", colnames(df))]        
treatment_columns <- colnames(df)[grepl("^treatment_", colnames(df))]        

# 删除第一行
df <- df[-1, ]
# 查看修改后的数据框
head(df)
write.csv(df, file = "./input_rename/36other_tsrna.csv")
```

3. 找到每个mirna的treatment_vs_control中位数，计算FC
```bash
cd /mnt/d/adult_dep/DE2
# 记得加“tsrna”
python median_tsrna.py ./input_rename/36other_tsrna.csv ./36other_tsrna/median.csv
```


② 缺失值替换为NA
```r
setwd("D:/adult_dep/DE2")
data <- read.csv("./input_rename/36other_tsrna.csv") 
data[data == 0] <- NA
write.csv(data, "./36other_tsrna/tsrna_NA.csv", row.names = FALSE)  
```



③ Wilcoxon秩和检验（相关样本）  

秩和检验（Wilcoxon秩和检验）是一种用于比较两组相关或无关样本的非参数检验方法。它不要求数据符合正态分布，因此适用于那些无法满足正态性假设的情况。Wilcoxon秩和检验（相关样本）： 适用于配对设计，比如同一组个体在两个不同时间点的观测。在这种情况下，对每对相关样本进行秩和的比较。    
* 原理：在Wilcoxon秩和检验中，基本思想是将每个组的观测值按大小排序，然后为每个观测值分配一个秩次，最后通过比较秩和来判断两组样本是否来自同一分布。这个检验的零假设是两组样本来自相同的分布。  


```r
# 设置工作目录
setwd("D:/adult_dep/DE2")

# 读取数据
data <- read.csv("./36other_tsrna/tsrna_NA.csv")

# 加载dplyr包
library(dplyr)

# 获取控制组和治疗组的列
control_columns <- grep("^control_", names(data), value = TRUE)
treatment_columns <- grep("^treatment_", names(data), value = TRUE)

# 执行 Wilcoxon 秩和检验
p_values_matrix <- matrix(NA, nrow = length(data$tsRNA), ncol = 3)
colnames(p_values_matrix) <- c("tsRNA", "p_value", "Q_value")

for (i in seq_along(data$tsRNA)) {
  tsRNA <- data$tsRNA[i]
  control_values <- unlist(data[data$tsRNA == tsRNA, control_columns])
  treatment_values <- unlist(data[data$tsRNA == tsRNA, treatment_columns])
  
  # 执行 Wilcoxon 秩和检验
  result <- try(wilcox.test(treatment_values, control_values))

  # 存储 p-value
  if (!inherits(result, "try-error")) {
    p_values_matrix[i, 1] <- tsRNA
    p_values_matrix[i, 2] <- result$p.value
  } else {
    cat("Error for tsRNA:", tsRNA, "\n")
  }
}

# 使用 Benjamini-Hochberg 方法进行多重比较校正
p_values_matrix[, 3] <- p.adjust(p_values_matrix[, 2], method = "BH")

# 输出结果
print(p_values_matrix)
# 将 p_values_matrix 保存到 CSV 文件
write.csv(p_values_matrix, file = "./36other_tsrna/de_values_matrix.csv", row.names = FALSE)
```

④ 合并文件
```r
library(dplyr)

median_data <- read.csv("./36other_tsrna/median.csv")
de_values_data <- read.csv("./36other_tsrna/de_values_matrix.csv")
# 基于 miRNA 列合并数据
merged_data <- inner_join(median_data, de_values_data, by = "tsRNA")

# 查看合并后的数据
print(merged_data)
write.csv(merged_data, file = "./36other_tsrna/result.csv", row.names = FALSE)

diff_gene <- subset(merged_data, Q_value < 0.05 & abs(log2fold_change) > 0.58)
write.csv(diff_gene, file="./36other_tsrna/sig_q0.05_fc15.csv", quote = F)
```



⑤ 绘图
```r
# 加载包
library(ggplot2)
library(pheatmap)

merged_data$lgpadj <- -log10(merged_data$Q_value)

# 绘制火山图
volcano_plot <- ggplot(merged_data, aes(x = log2fold_change, y = lgpadj)) +
  geom_point(aes(color = ifelse(Q_value < 0.05 & log2fold_change > 1, "Up", 
                                ifelse(Q_value < 0.05 & log2fold_change < -1, "Down", "Not Sig"))), size = 2) +
  scale_color_manual(values = c("Up" = "#D32F2F", "Down" = "#1976D2", "Not Sig" = "grey70")) +
  labs(x = "Log2 Fold Change", y = "-Log10(padj)", title = "Volcano Plot") +
  theme_minimal(base_size = 14) +
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(),  # 去掉主网格线
        panel.grid.minor = element_blank(),  # 去掉次网格线
        axis.line = element_line(color = "black")) + # 显示横纵坐标轴
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", size = 0.5) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black", size = 0.5)

ggsave("./bp_tsrna/volcano_plot.pdf", plot = volcano_plot, device = "pdf", dpi = 600, width = 8, height = 6)
```

```r
library(ggplot2)
library(pheatmap)

merged_data <- read.csv("./36other_tsrna/result.csv")
merged_data$lgpadj <- -log10(merged_data$Q_value)

# Define categories for points
merged_data$category <- ifelse(merged_data$Q_value < 0.05 & merged_data$log2fold_change > 0.58, "Up", 
                                ifelse(merged_data$Q_value < 0.05 & merged_data$log2fold_change < -0.58, "Down", "Not Sig"))

# Calculate counts for each category
category_counts <- table(merged_data$category)

# Create a volcano plot
volcano_plot <- ggplot(merged_data, aes(x = log2fold_change, y = lgpadj)) +
  geom_point(aes(color = category), size = 2) +
  scale_color_manual(values = c("Up" = "#D32F2F", "Down" = "#1976D2", "Not Sig" = "grey70")) +
  labs(x = "Log2 Fold Change", y = "-Log10(padj)", title = "Volcano Plot") +
  theme_minimal(base_size = 14) +
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        axis.line = element_line(color = "black"),  # Show axis lines
        axis.ticks = element_line(color = "black"),  # Show axis ticks
        axis.ticks.length = unit(0.25, "cm")) +  # Set tick length
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", size = 0.5) +
  geom_vline(xintercept = c(-0.58, 0.58), linetype = "dashed", color = "black", size = 0.5) +
  
  # Add category count labels
  geom_text(aes(x = 1, y = 30, label = paste("Up:", category_counts["Up"])), color = "#D32F2F", size = 5, hjust = 0) +
  geom_text(aes(x = 1, y = 28, label = paste("Down:", category_counts["Down"])), color = "#1976D2", size = 5, hjust = 0) +
  geom_text(aes(x = 1, y = 26, label = paste("Not Sig:", category_counts["Not Sig"])), color = "grey70", size = 5, hjust = 0)

# Save the plot
ggsave("./36other_tsrna/36other_tsrna_volcano_plot.pdf", plot = volcano_plot, device = "pdf", dpi = 600, width = 8, height = 6)
```




* 图中标注

```r
# 图中标注名称
annotations <- data.frame(
  tsRNA = c("ACTTGACCGCTCTGACCA", "ATCCTGCCGACTACGCCA", "GCATTGGTGGTTCAGTGGTAGAATTCT","GCCCGGATAGCTCAGTCGGTAGAG","GCGGGAGACCGGGGTTCGATTCCCCGACGGGGAGC","GTTTCCGTAGTGTAGTGGTTATC"),
  new_label = c("tsRNA-#9", "tsRNA-#17", "tsRNA-#36","tsRNA-#45","tsRNA-#51","tsRNA-#56")
)

# 提取需要标注的 miRNA 数据
annotations_data <- merged_data[merged_data$tsRNA %in% annotations$tsRNA, ]

# 创建火山图
volcano_plot <- ggplot(merged_data, aes(x = log2fold_change, y = lgpadj)) +
  geom_point(aes(color = ifelse(Q_value < 0.05 & log2fold_change > 0.58, "Up", 
                                ifelse(Q_value < 0.05 & log2fold_change < -0.58, "Down", "Not Sig"))), size = 2) +
  scale_color_manual(values = c("Up" = "#D32F2F", "Down" = "#1976D2", "Not Sig" = "grey70")) +
  labs(x = "Log2 Fold Change", y = "-Log10(padj)", title = "Volcano Plot") +
  theme_minimal(base_size = 14) +
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", size = 0.5) +
  geom_vline(xintercept = c(-0.58, 0.58), linetype = "dashed", color = "black", size = 0.5) +
  geom_point(data = annotations_data, aes(color = "Highlighted"), size = 3, shape = 1, stroke = 1.5, fill = NA) +  # 添加圈圈
  geom_text(data = merge(annotations_data, annotations, by = "tsRNA"), aes(label = new_label), 
            vjust = -1, hjust = 0.5, size = 3, color = "black")  # 使用新的标签


ggsave("./bp_tsrna/label_bp_tsrna_volcano_plot.pdf", plot = volcano_plot, device = "pdf", dpi = 600, width = 8, height = 6)

```

```r
# 绘制热图
library(ggplot2)
library(tidyr)
library(dplyr)

# 读取数据
merged_data <- read.csv("./bp_tsrna/result.csv")

# 指定需要标注的 tsRNA 和对应的新名称
tsRNAs_to_annotate <- c("ACTTGACCGCTCTGACCA", "ATCCTGCCGACTACGCCA", "GCATTGGTGGTTCAGTGGTAGAATTCT","GCCCGGATAGCTCAGTCGGTAGAG","GCGGGAGACCGGGGTTCGATTCCCCGACGGGGAGC","GTTTCCGTAGTGTAGTGGTTATC")

new_names <- c("tsRNA-#9", "tsRNA-#17", "tsRNA-#36","tsRNA-#45","tsRNA-#51","tsRNA-#56")

# 创建一个命名映射
name_mapping <- setNames(new_names, tsRNAs_to_annotate)

# 提取需要的 tsRNA 数据并设置顺序
heatmap_data <- merged_data %>%
  filter(tsRNA %in% tsRNAs_to_annotate) %>%
  select(tsRNA, log2control_median, log2treatment_median) %>%
  pivot_longer(cols = starts_with("log2"), names_to = "Condition", values_to = "Expression") %>%
  mutate(tsRNA = factor(name_mapping[tsRNA], levels = new_names))  # 使用新名称

# 定义颜色从白色到深蓝色
start_color <- rgb(255/255, 255/255, 255/255)  # 白色
end_color <- rgb(0/255, 51/255, 102/255)        # 深蓝色

# 生成渐变
num_colors <- 10
gradient_colors <- colorRampPalette(c(start_color, end_color))(num_colors)

# 绘制热图
heatmap_plot <- ggplot(heatmap_data, aes(x = tsRNA, y = Condition, fill = Expression)) +
  geom_tile(color = "white") +
  scale_fill_gradientn(colors = gradient_colors, 
                       name = "Expression Level") +
  labs(x = "tsRNA", y = "Condition", title = "Heatmap of tsRNA Expression") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 保存热图
ggsave("./bp_tsrna/5tsRNA_heatmap.pdf", plot = heatmap_plot, device = "pdf", dpi = 600, width = 8, height = 6)
```










































## BD vs HC

### 1. bp_mirna

2. 标注control_vs_treatment，手动标注
```r
setwd("D:/adult_dep/DE2")
df <- read.csv("./input/rpm_50%_log_mirna.csv", header = TRUE, row.names = 1)
df <- df[, df[1, ] %in% c(0, 2)]
colnames(df) <- ifelse(df[1, ] == 0, 
                       paste0("control_", colnames(df)), 
                       paste0("treatment_", colnames(df)))
control_columns <- colnames(df)[grepl("^control_", colnames(df))]        

# 删除第一行
df <- df[-1, ]
# 查看修改后的数据框
head(df)
write.csv(df, file = "./input_rename/bp_hc_mirna.csv")
```

3. 找到每个mirna的treatment_vs_control中位数，计算FC
```bash
cd /mnt/d/adult_dep/DE2
# 记得加“mirna”
python median_mirna.py ./input_rename/bp_hc_mirna.csv ./bp_hc_mirna/median.csv

```


② 缺失值替换为NA
```r
setwd("D:/adult_dep/DE2")
data <- read.csv("./input_rename/bp_hc_mirna.csv") 
data[data == 0] <- NA
write.csv(data, "./bp_hc_mirna/mirna_NA.csv", row.names = FALSE)  
```



③ Wilcoxon秩和检验（相关样本）  

秩和检验（Wilcoxon秩和检验）是一种用于比较两组相关或无关样本的非参数检验方法。它不要求数据符合正态分布，因此适用于那些无法满足正态性假设的情况。Wilcoxon秩和检验（相关样本）： 适用于配对设计，比如同一组个体在两个不同时间点的观测。在这种情况下，对每对相关样本进行秩和的比较。    
* 原理：在Wilcoxon秩和检验中，基本思想是将每个组的观测值按大小排序，然后为每个观测值分配一个秩次，最后通过比较秩和来判断两组样本是否来自同一分布。这个检验的零假设是两组样本来自相同的分布。  


```r
# 设置工作目录
setwd("D:/adult_dep/DE2")

# 读取数据
data <- read.csv("./bp_hc_mirna/mirna_NA.csv")

# 加载dplyr包
library(dplyr)

# 获取控制组和治疗组的列
control_columns <- grep("^control_", names(data), value = TRUE)
treatment_columns <- grep("^treatment_", names(data), value = TRUE)


# 执行 Wilcoxon 秩和检验
p_values_matrix <- matrix(NA, nrow = length(data$miRNA), ncol = 3)
colnames(p_values_matrix) <- c("miRNA", "p_value", "Q_value")

for (i in seq_along(data$miRNA)) {
  miRNA <- data$miRNA[i]
  control_values <- unlist(data[data$miRNA == miRNA, control_columns])
  treatment_values <- unlist(data[data$miRNA == miRNA, treatment_columns])
  
  # 执行 Wilcoxon 秩和检验
  result <- try(wilcox.test(treatment_values, control_values))

  # 存储 p-value
  if (!inherits(result, "try-error")) {
    p_values_matrix[i, 1] <- miRNA
    p_values_matrix[i, 2] <- result$p.value
  } else {
    cat("Error for miRNA:", miRNA, "\n")
  }
}

# 使用 Benjamini-Hochberg 方法进行多重比较校正
p_values_matrix[, 3] <- p.adjust(p_values_matrix[, 2], method = "BH")

# 输出结果
print(p_values_matrix)
# 将 p_values_matrix 保存到 CSV 文件
write.csv(p_values_matrix, file = "./bp_hc_mirna/de_values_matrix.csv", row.names = FALSE)
```

④ 合并文件
```r
library(dplyr)

median_data <- read.csv("./bp_hc_mirna/median.csv")
de_values_data <- read.csv("./bp_hc_mirna/de_values_matrix.csv")
# 基于 miRNA 列合并数据
merged_data <- inner_join(median_data, de_values_data, by = "miRNA")

# 查看合并后的数据
print(merged_data)
write.csv(merged_data, file = "./bp_hc_mirna/result.csv", row.names = FALSE)
```

```r
library(dplyr)

median_data <- read.csv("./bp_hc_mirna/median.csv")
de_values_data <- read.csv("./bp_hc_mirna/de_values_matrix.csv")
# 基于 miRNA 列合并数据
merged_data <- inner_join(median_data, de_values_data, by = "miRNA")


diff_gene <- subset(merged_data, Q_value < 0.05 & abs(log2fold_change) > 0.58)
write.csv(diff_gene, file="./bp_hc_mirna/sig_q0.05_fc15.csv", quote = F)
```









⑤ 绘图

```r
# 加载包
library(ggplot2)
library(pheatmap)
setwd("D:/adult_dep/DE2")

merged_data <- read.csv("./bp_hc_mirna/result.csv")
merged_data$lgpadj <- -log10(merged_data$Q_value)

# Define categories for points
merged_data$category <- ifelse(merged_data$Q_value < 0.05 & merged_data$log2fold_change > 0.58, "Up", 
                                ifelse(merged_data$Q_value < 0.05 & merged_data$log2fold_change < -0.58, "Down", "Not Sig"))

# Calculate counts for each category
category_counts <- table(merged_data$category)

# Create a volcano plot
volcano_plot <- ggplot(merged_data, aes(x = log2fold_change, y = lgpadj)) +
  geom_point(aes(color = category), size = 2) +
  scale_color_manual(values = c("Up" = "#D32F2F", "Down" = "#1976D2", "Not Sig" = "grey70")) +
  labs(x = "Log2 Fold Change", y = "-Log10(padj)", title = "Volcano Plot") +
  theme_minimal(base_size = 14) +
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        axis.line = element_line(color = "black"),  # Show axis lines
        axis.ticks = element_line(color = "black"),  # Show axis ticks
        axis.ticks.length = unit(0.25, "cm")) +  # Set tick length
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", size = 0.5) +
  geom_vline(xintercept = c(-0.58, 0.58), linetype = "dashed", color = "black", size = 0.5) +
  
  # Add category count labels
  geom_text(aes(x = 2, y = 30, label = paste("Up:", category_counts["Up"])), color = "#D32F2F", size = 5, hjust = 0) +
  geom_text(aes(x = 2, y = 28, label = paste("Down:", category_counts["Down"])), color = "#1976D2", size = 5, hjust = 0) +
  geom_text(aes(x = 2, y = 26, label = paste("Not Sig:", category_counts["Not Sig"])), color = "grey70", size = 5, hjust = 0)

ggsave("./bp_hc_mirna/bp_hc_mirna_volcano_plot.pdf", plot = volcano_plot, device = "pdf", dpi = 600, width = 8, height = 6)
```

```r
# 图中标注名称
miRNAs_to_annotate <- c("hsa-mir-28_hsa-miR-28-3p", "hsa-mir-410_hsa-miR-410-3p", 
                         "hsa-mir-4433b_hsa-miR-4433b-5p")

# 提取需要标注的 miRNA 数据
annotations <- merged_data[merged_data$miRNA %in% miRNAs_to_annotate, ]

# 创建火山图
volcano_plot <- ggplot(merged_data, aes(x = log2fold_change, y = lgpadj)) +
  geom_point(aes(color = ifelse(Q_value < 0.05 & log2fold_change > 0.58, "Up", 
                                ifelse(Q_value < 0.05 & log2fold_change < -0.58, "Down", "Not Sig"))), size = 2) +
  scale_color_manual(values = c("Up" = "#D32F2F", "Down" = "#1976D2", "Not Sig" = "grey70")) +
  labs(x = "Log2 Fold Change", y = "-Log10(padj)", title = "Volcano Plot") +
  theme_minimal(base_size = 14) +
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", size = 0.5) +
  geom_vline(xintercept = c(-0.58, 0.58), linetype = "dashed", color = "black", size = 0.5) +
  geom_point(data = annotations, aes(color = "Highlighted"), size = 3, shape = 1, stroke = 1.5, fill = NA) +  # 添加圈圈
  geom_text(data = annotations, aes(label = miRNA), vjust = -1, hjust = 0.5, size = 3, color = "black")


ggsave("./bp_hc_mirna/label_bp_hc_mirna_volcano_plot.pdf", plot = volcano_plot, device = "pdf", dpi = 600, width = 8, height = 6)

```











### 2. bp_hc_pirna

2. 标注control_vs_treatment，手动标注
```r
setwd("D:/adult_dep/DE2")
df <- read.csv("./input/rpm_50%_log_pirna.csv", header = TRUE, row.names = 1)
df <- df[, df[1, ] %in% c(0, 2)]
colnames(df) <- ifelse(df[1, ] == 0, 
                       paste0("control_", colnames(df)), 
                       paste0("treatment_", colnames(df)))
control_columns <- colnames(df)[grepl("^control_", colnames(df))]        

# 删除第一行
df <- df[-1, ]
# 查看修改后的数据框
head(df)
write.csv(df, file = "./input_rename/bp_hc_pirna.csv")
```

3. 找到每个mirna的treatment_vs_control中位数，计算FC
```bash
cd /mnt/d/adult_dep/DE2
# 记得加“pirna”
python median_pirna.py ./input_rename/bp_hc_pirna.csv ./bp_hc_pirna/median.csv
```


② 缺失值替换为NA
```r
setwd("D:/adult_dep/DE2")
data <- read.csv("./input_rename/bp_hc_pirna.csv") 
data[data == 0] <- NA
write.csv(data, "./bp_hc_pirna/pirna_NA.csv", row.names = FALSE)  
```



③ Wilcoxon秩和检验（相关样本）  

秩和检验（Wilcoxon秩和检验）是一种用于比较两组相关或无关样本的非参数检验方法。它不要求数据符合正态分布，因此适用于那些无法满足正态性假设的情况。Wilcoxon秩和检验（相关样本）： 适用于配对设计，比如同一组个体在两个不同时间点的观测。在这种情况下，对每对相关样本进行秩和的比较。    
* 原理：在Wilcoxon秩和检验中，基本思想是将每个组的观测值按大小排序，然后为每个观测值分配一个秩次，最后通过比较秩和来判断两组样本是否来自同一分布。这个检验的零假设是两组样本来自相同的分布。  


```r
# 设置工作目录
setwd("D:/adult_dep/DE2")

# 读取数据
data <- read.csv("./bp_hc_pirna/pirna_NA.csv")

# 加载dplyr包
library(dplyr)

# 获取控制组和治疗组的列
control_columns <- grep("^control_", names(data), value = TRUE)
treatment_columns <- grep("^treatment_", names(data), value = TRUE)
# 835

# 执行 Wilcoxon 秩和检验
p_values_matrix <- matrix(NA, nrow = length(data$piRNA), ncol = 3)
colnames(p_values_matrix) <- c("piRNA", "p_value", "Q_value")

for (i in seq_along(data$piRNA)) {
  piRNA <- data$piRNA[i]
  control_values <- unlist(data[data$piRNA == piRNA, control_columns])
  treatment_values <- unlist(data[data$piRNA == piRNA, treatment_columns])
  
  # 执行 Wilcoxon 秩和检验
  result <- try(wilcox.test(treatment_values, control_values))

  # 存储 p-value
  if (!inherits(result, "try-error")) {
    p_values_matrix[i, 1] <- piRNA
    p_values_matrix[i, 2] <- result$p.value
  } else {
    cat("Error for piRNA:", piRNA, "\n")
  }
}

# 使用 Benjamini-Hochberg 方法进行多重比较校正
p_values_matrix[, 3] <- p.adjust(p_values_matrix[, 2], method = "BH")

# 输出结果
print(p_values_matrix)
# 将 p_values_matrix 保存到 CSV 文件
write.csv(p_values_matrix, file = "./bp_hc_pirna/de_values_matrix.csv", row.names = FALSE)
```

④ 合并文件
```r
library(dplyr)

median_data <- read.csv("./bp_hc_pirna/median.csv")
de_values_data <- read.csv("./bp_hc_pirna/de_values_matrix.csv")
# 基于 miRNA 列合并数据
merged_data <- inner_join(median_data, de_values_data, by = "piRNA")

# 查看合并后的数据
print(merged_data)
write.csv(merged_data, file = "./bp_hc_pirna/result.csv", row.names = FALSE)
```

```r
library(dplyr)
setwd("D:/adult_dep/DE2")

merged_data <- read.csv("./bp_hc_pirna/result.csv")

diff_gene <- subset(merged_data, Q_value < 0.05 & abs(log2fold_change) > 0.58)
write.csv(diff_gene, file="./bp_hc_pirna/sig_q0.05_fc15.csv", quote = F)
```




⑤ 绘图
```r
# 加载包
library(ggplot2)
library(pheatmap)

merged_data$lgpadj <- -log10(merged_data$Q_value)


# Define categories for points
merged_data$category <- ifelse(merged_data$Q_value < 0.05 & merged_data$log2fold_change > 0.58, "Up", 
                                ifelse(merged_data$Q_value < 0.05 & merged_data$log2fold_change < -0.58, "Down", "Not Sig"))

# Calculate counts for each category
category_counts <- table(merged_data$category)

# Create a volcano plot
volcano_plot <- ggplot(merged_data, aes(x = log2fold_change, y = lgpadj)) +
  geom_point(aes(color = category), size = 2) +
  scale_color_manual(values = c("Up" = "#D32F2F", "Down" = "#1976D2", "Not Sig" = "grey70")) +
  labs(x = "Log2 Fold Change", y = "-Log10(padj)", title = "Volcano Plot") +
  theme_minimal(base_size = 14) +
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        axis.line = element_line(color = "black"),  # Show axis lines
        axis.ticks = element_line(color = "black"),  # Show axis ticks
        axis.ticks.length = unit(0.25, "cm")) +  # Set tick length
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", size = 0.5) +
  geom_vline(xintercept = c(-0.58, 0.58), linetype = "dashed", color = "black", size = 0.5) +
  
  # Add category count labels
  geom_text(aes(x = 2, y = 30, label = paste("Up:", category_counts["Up"])), color = "#D32F2F", size = 5, hjust = 0) +
  geom_text(aes(x = 2, y = 28, label = paste("Down:", category_counts["Down"])), color = "#1976D2", size = 5, hjust = 0) +
  geom_text(aes(x = 2, y = 26, label = paste("Not Sig:", category_counts["Not Sig"])), color = "grey70", size = 5, hjust = 0)

ggsave("./bp_hc_pirna/volcano_plot.pdf", plot = volcano_plot, device = "pdf", dpi = 600, width = 8, height = 6)
```




* 图中标注

```r
# 图中标注名称
miRNAs_to_annotate <- c("piR-hsa-135526", "piR-hsa-147698","piR-hsa-151782","piR-hsa-184607","piR-hsa-789156")

# 提取需要标注的 miRNA 数据
annotations <- merged_data[merged_data$piRNA %in% miRNAs_to_annotate, ]

# 创建火山图
volcano_plot <- ggplot(merged_data, aes(x = log2fold_change, y = lgpadj)) +
  geom_point(aes(color = ifelse(Q_value < 0.05 & log2fold_change > 0.58, "Up", 
                                ifelse(Q_value < 0.05 & log2fold_change < -0.58, "Down", "Not Sig"))), size = 2) +
  scale_color_manual(values = c("Up" = "#D32F2F", "Down" = "#1976D2", "Not Sig" = "grey70")) +
  labs(x = "Log2 Fold Change", y = "-Log10(padj)", title = "Volcano Plot") +
  theme_minimal(base_size = 14) +
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", size = 0.5) +
  geom_vline(xintercept = c(-0.58, 0.58), linetype = "dashed", color = "black", size = 0.5) +
  geom_point(data = annotations, aes(color = "Highlighted"), size = 3, shape = 1, stroke = 1.5, fill = NA) +  # 添加圈圈
  geom_text(data = annotations, aes(label = piRNA), vjust = -1, hjust = 0.5, size = 3, color = "black")


ggsave("./bp_hc_pirna/label_bp_hc_pirna_volcano_plot.pdf", plot = volcano_plot, device = "pdf", dpi = 600, width = 8, height = 6)

```

```r
# 绘制热图
library(ggplot2)
library(tidyr)
library(dplyr)

# 读取数据
merged_data <- read.csv("./bp_hc_pirna/result.csv")

# 指定要绘制热图的 miRNA
piRNAs_to_annotate <- c("piR-hsa-135526", "piR-hsa-147698","piR-hsa-151782","piR-hsa-184607","piR-hsa-789156")

# 提取需要的 piRNA 数据并设置顺序
heatmap_data <- merged_data %>%
  filter(piRNA %in% piRNAs_to_annotate) %>%
  select(piRNA, log2control_median, log2treatment_median) %>%
  pivot_longer(cols = starts_with("log2"), names_to = "Condition", values_to = "Expression") %>%
  mutate(piRNA = factor(piRNA, levels = piRNAs_to_annotate))  # 设置因子顺序

# 定义颜色从白色到深蓝色
start_color <- rgb(255/255, 255/255, 255/255)  # 白色
end_color <- rgb(0/255, 51/255, 102/255)        # 深蓝色

# 生成渐变
num_colors <- 10
gradient_colors <- colorRampPalette(c(start_color, end_color))(num_colors)

# 绘制热图
heatmap_plot <- ggplot(heatmap_data, aes(x = piRNA, y = Condition, fill = Expression)) +
  geom_tile(color = "white") +
  scale_fill_gradientn(colors = gradient_colors, 
                       name = "Expression Level") +
  labs(x = "piRNA", y = "Condition", title = "Heatmap of piRNA Expression") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 保存热图
ggsave("./bp_hc_pirna/5piRNA_heatmap.pdf", plot = heatmap_plot, device = "pdf", dpi = 600, width = 8, height = 6)
```






## 3. bp_hc_rsrna

2. 标注control_vs_treatment，手动标注
```r
setwd("D:/adult_dep/DE2")
df <- read.csv("./input/rpm_50%_log_rsrna.csv", header = TRUE, row.names = 1)
df <- df[, df[1, ] %in% c(0, 2)]
colnames(df) <- ifelse(df[1, ] == 0, 
                       paste0("control_", colnames(df)), 
                       paste0("treatment_", colnames(df)))
control_columns <- colnames(df)[grepl("^control_", colnames(df))]       

# 删除第一行
df <- df[-1, ]
# 查看修改后的数据框
head(df)
write.csv(df, file = "./input_rename/bp_hc_rsrna.csv")
```

3. 找到每个mirna的treatment_vs_control中位数，计算FC
```bash
cd /mnt/d/adult_dep/DE2
# 记得加“rsrna”
python median_rsrna.py ./input_rename/bp_hc_rsrna.csv ./bp_hc_rsrna/median.csv
```


② 缺失值替换为NA
```r
setwd("D:/adult_dep/DE2")
data <- read.csv("./input_rename/bp_hc_rsrna.csv") 
data[data == 0] <- NA
write.csv(data, "./bp_hc_rsrna/rsrna_NA.csv", row.names = FALSE)  
```



③ Wilcoxon秩和检验（相关样本）  

秩和检验（Wilcoxon秩和检验）是一种用于比较两组相关或无关样本的非参数检验方法。它不要求数据符合正态分布，因此适用于那些无法满足正态性假设的情况。Wilcoxon秩和检验（相关样本）： 适用于配对设计，比如同一组个体在两个不同时间点的观测。在这种情况下，对每对相关样本进行秩和的比较。    
* 原理：在Wilcoxon秩和检验中，基本思想是将每个组的观测值按大小排序，然后为每个观测值分配一个秩次，最后通过比较秩和来判断两组样本是否来自同一分布。这个检验的零假设是两组样本来自相同的分布。  


```r
# 设置工作目录
setwd("D:/adult_dep/DE2")

# 读取数据
data <- read.csv("./bp_hc_rsrna/rsrna_NA.csv")

# 加载dplyr包
library(dplyr)

# 获取控制组和治疗组的列
control_columns <- grep("^control_", names(data), value = TRUE)
treatment_columns <- grep("^treatment_", names(data), value = TRUE)

# 执行 Wilcoxon 秩和检验
p_values_matrix <- matrix(NA, nrow = length(data$rsRNA), ncol = 3)
colnames(p_values_matrix) <- c("rsRNA", "p_value", "Q_value")

for (i in seq_along(data$rsRNA)) {
  rsRNA <- data$rsRNA[i]
  control_values <- unlist(data[data$rsRNA == rsRNA, control_columns])
  treatment_values <- unlist(data[data$rsRNA == rsRNA, treatment_columns])
  
  # 执行 Wilcoxon 秩和检验
  result <- try(wilcox.test(treatment_values, control_values))

  # 存储 p-value
  if (!inherits(result, "try-error")) {
    p_values_matrix[i, 1] <- rsRNA
    p_values_matrix[i, 2] <- result$p.value
  } else {
    cat("Error for rsRNA:", rsRNA, "\n")
  }
}

# 使用 Benjamini-Hochberg 方法进行多重比较校正
p_values_matrix[, 3] <- p.adjust(p_values_matrix[, 2], method = "BH")

# 输出结果
print(p_values_matrix)
write.csv(p_values_matrix, file = "./bp_hc_rsrna/de_values_matrix.csv", row.names = FALSE)
```

④ 合并文件
```r
library(dplyr)

median_data <- read.csv("./bp_hc_rsrna/median.csv")
de_values_data <- read.csv("./bp_hc_rsrna/de_values_matrix.csv")
# 基于 miRNA 列合并数据
merged_data <- inner_join(median_data, de_values_data, by = "rsRNA")

# 查看合并后的数据
print(merged_data)
write.csv(merged_data, file = "./bp_hc_rsrna/result.csv", row.names = FALSE)

library(dplyr)
setwd("D:/adult_dep/DE2")

merged_data <- read.csv("./bp_hc_rsrna/result.csv")

diff_gene <- subset(merged_data, Q_value < 0.05 & abs(log2fold_change) > 0.58)
write.csv(diff_gene, file="./bp_hc_rsrna/sig_q0.05_fc15.csv", quote = F)
```



⑤ 绘图

```r
# 加载包
library(ggplot2)
library(pheatmap)

merged_data <- read.csv("./bp_hc_rsrna/result.csv")
merged_data$lgpadj <- -log10(merged_data$Q_value)

# Define categories for points
merged_data$category <- ifelse(merged_data$Q_value < 0.05 & merged_data$log2fold_change > 0.58, "Up", 
                                ifelse(merged_data$Q_value < 0.05 & merged_data$log2fold_change < -0.58, "Down", "Not Sig"))

# Calculate counts for each category
category_counts <- table(merged_data$category)

# Create a volcano plot
volcano_plot <- ggplot(merged_data, aes(x = log2fold_change, y = lgpadj)) +
  geom_point(aes(color = category), size = 2) +
  scale_color_manual(values = c("Up" = "#D32F2F", "Down" = "#1976D2", "Not Sig" = "grey70")) +
  labs(x = "Log2 Fold Change", y = "-Log10(padj)", title = "Volcano Plot") +
  theme_minimal(base_size = 14) +
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        axis.line = element_line(color = "black"),  # Show axis lines
        axis.ticks = element_line(color = "black"),  # Show axis ticks
        axis.ticks.length = unit(0.25, "cm")) +  # Set tick length
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", size = 0.5) +
  geom_vline(xintercept = c(-0.58, 0.58), linetype = "dashed", color = "black", size = 0.5) +
  
  # Add category count labels
  geom_text(aes(x = 2, y = 30, label = paste("Up:", category_counts["Up"])), color = "#D32F2F", size = 5, hjust = 0) +
  geom_text(aes(x = 2, y = 28, label = paste("Down:", category_counts["Down"])), color = "#1976D2", size = 5, hjust = 0) +
  geom_text(aes(x = 2, y = 26, label = paste("Not Sig:", category_counts["Not Sig"])), color = "grey70", size = 5, hjust = 0)

ggsave("./bp_hc_rsrna/bp_hc_rsrna_volcano_plot.pdf", plot = volcano_plot, device = "pdf", dpi = 600, width = 8, height = 6)
```




* 图中标注

```r
# 图中标注名称
annotations <- data.frame(
  rsRNA = c("ACAATGTAGGTAAGGGAAGTCG", "ACGGACCAAGGAGTCTAACACGTGCGCG", "CGACTCTTAGCGGTGGATCACTCGGCTCGTGCG","GGAACAATGTAGGTAAGGGA","GGAGTCTAACACGTGCGCG","TGGGAGACCGCCTGGGAATACCGGGTGCTGTAGGCTT", "GACTCTTAGCGGTGGATCACTCGGCTCGTG"),
  new_label = c("rsRNA-#18", "rsRNA-#29", "rsRNA-#95","rsRNA-#176","rsRNA-#182","rsRNA-#261", "rsRNA-#155")
)

# 提取需要标注的 miRNA 数据
annotations_data <- merged_data[merged_data$rsRNA %in% annotations$rsRNA, ]

# 创建火山图
volcano_plot <- ggplot(merged_data, aes(x = log2fold_change, y = lgpadj)) +
  geom_point(aes(color = ifelse(Q_value < 0.05 & log2fold_change > 0.58, "Up", 
                                ifelse(Q_value < 0.05 & log2fold_change < -0.58, "Down", "Not Sig"))), size = 2) +
  scale_color_manual(values = c("Up" = "#D32F2F", "Down" = "#1976D2", "Not Sig" = "grey70")) +
  labs(x = "Log2 Fold Change", y = "-Log10(padj)", title = "Volcano Plot") +
  theme_minimal(base_size = 14) +
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", size = 0.5) +
  geom_vline(xintercept = c(-0.58, 0.58), linetype = "dashed", color = "black", size = 0.5) +
  geom_point(data = annotations_data, aes(color = "Highlighted"), size = 3, shape = 1, stroke = 1.5, fill = NA) +  # 添加圈圈
  geom_text(data = merge(annotations_data, annotations, by = "rsRNA"), aes(label = new_label), 
            vjust = -1, hjust = 0.5, size = 3, color = "black")  # 使用新的标签


ggsave("./bp_hc_rsrna/label_bp_hc_rsrna_volcano_plot.pdf", plot = volcano_plot, device = "pdf", dpi = 600, width = 8, height = 6)

```

```r
# 绘制热图
library(ggplot2)
library(tidyr)
library(dplyr)

# 读取数据
merged_data <- read.csv("./bp_hc_rsrna/result.csv")

# 指定需要标注的 tsRNA 和对应的新名称
rsRNAs_to_annotate <- c("ACAATGTAGGTAAGGGAAGTCG", "ACGGACCAAGGAGTCTAACACGTGCGCG", "CGACTCTTAGCGGTGGATCACTCGGCTCGTGCG","GGAACAATGTAGGTAAGGGA","GGAGTCTAACACGTGCGCG","TGGGAGACCGCCTGGGAATACCGGGTGCTGTAGGCTT", "GACTCTTAGCGGTGGATCACTCGGCTCGTG")

new_names <- c("rsRNA-#18", "rsRNA-#29", "rsRNA-#95","rsRNA-#176","rsRNA-#182","rsRNA-#261", "rsRNA-#155")

# 创建一个命名映射
name_mapping <- setNames(new_names, rsRNAs_to_annotate)

# 提取需要的 rsRNA 数据并设置顺序
heatmap_data <- merged_data %>%
  filter(rsRNA %in% rsRNAs_to_annotate) %>%
  select(rsRNA, log2control_median, log2treatment_median) %>%
  pivot_longer(cols = starts_with("log2"), names_to = "Condition", values_to = "Expression") %>%
  mutate(rsRNA = factor(name_mapping[rsRNA], levels = new_names))  # 使用新名称

# 定义颜色从白色到深蓝色
start_color <- rgb(255/255, 255/255, 255/255)  # 白色
end_color <- rgb(0/255, 51/255, 102/255)        # 深蓝色

# 生成渐变
num_colors <- 10
gradient_colors <- colorRampPalette(c(start_color, end_color))(num_colors)

# 绘制热图
heatmap_plot <- ggplot(heatmap_data, aes(x = rsRNA, y = Condition, fill = Expression)) +
  geom_tile(color = "white") +
  scale_fill_gradientn(colors = gradient_colors, 
                       name = "Expression Level") +
  labs(x = "rsRNA", y = "Condition", title = "Heatmap of rsRNA Expression") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 保存热图
ggsave("./bp_hc_rsrna/5bp_hc_heatmap.pdf", plot = heatmap_plot, device = "pdf", dpi = 600, width = 8, height = 6)
```





## 4. bp_hc_tsrna

2. 标注control_vs_treatment，手动标注
```r
setwd("D:/adult_dep/DE2")
df <- read.csv("./input/rpm_50%_log_tsrna.csv", header = TRUE, row.names = 1)
df <- df[, df[1, ] %in% c(0, 2)]
colnames(df) <- ifelse(df[1, ] == 0, 
                       paste0("control_", colnames(df)), 
                       paste0("treatment_", colnames(df)))
control_columns <- colnames(df)[grepl("^control_", colnames(df))]          

# 删除第一行
df <- df[-1, ]
# 查看修改后的数据框
head(df)
write.csv(df, file = "./input_rename/bp_hc_tsrna.csv")
```

3. 找到每个mirna的treatment_vs_control中位数，计算FC
```bash
cd /mnt/d/adult_dep/DE2
# 记得加“tsrna”
python median_tsrna.py ./input_rename/bp_hc_tsrna.csv ./bp_hc_tsrna/median.csv
```


② 缺失值替换为NA
```r
setwd("D:/adult_dep/DE2")
data <- read.csv("./input_rename/bp_hc_tsrna.csv") 
data[data == 0] <- NA
write.csv(data, "./bp_hc_tsrna/tsrna_NA.csv", row.names = FALSE)  
```



③ Wilcoxon秩和检验（相关样本）  

秩和检验（Wilcoxon秩和检验）是一种用于比较两组相关或无关样本的非参数检验方法。它不要求数据符合正态分布，因此适用于那些无法满足正态性假设的情况。Wilcoxon秩和检验（相关样本）： 适用于配对设计，比如同一组个体在两个不同时间点的观测。在这种情况下，对每对相关样本进行秩和的比较。    
* 原理：在Wilcoxon秩和检验中，基本思想是将每个组的观测值按大小排序，然后为每个观测值分配一个秩次，最后通过比较秩和来判断两组样本是否来自同一分布。这个检验的零假设是两组样本来自相同的分布。  


```r
# 设置工作目录
setwd("D:/adult_dep/DE2")

# 读取数据
data <- read.csv("./bp_hc_tsrna/tsrna_NA.csv")

# 加载dplyr包
library(dplyr)

# 获取控制组和治疗组的列
control_columns <- grep("^control_", names(data), value = TRUE)
treatment_columns <- grep("^treatment_", names(data), value = TRUE)

# 执行 Wilcoxon 秩和检验
p_values_matrix <- matrix(NA, nrow = length(data$tsRNA), ncol = 3)
colnames(p_values_matrix) <- c("tsRNA", "p_value", "Q_value")

for (i in seq_along(data$tsRNA)) {
  tsRNA <- data$tsRNA[i]
  control_values <- unlist(data[data$tsRNA == tsRNA, control_columns])
  treatment_values <- unlist(data[data$tsRNA == tsRNA, treatment_columns])
  
  # 执行 Wilcoxon 秩和检验
  result <- try(wilcox.test(treatment_values, control_values))

  # 存储 p-value
  if (!inherits(result, "try-error")) {
    p_values_matrix[i, 1] <- tsRNA
    p_values_matrix[i, 2] <- result$p.value
  } else {
    cat("Error for tsRNA:", tsRNA, "\n")
  }
}

# 使用 Benjamini-Hochberg 方法进行多重比较校正
p_values_matrix[, 3] <- p.adjust(p_values_matrix[, 2], method = "BH")

# 输出结果
print(p_values_matrix)
# 将 p_values_matrix 保存到 CSV 文件
write.csv(p_values_matrix, file = "./bp_hc_tsrna/de_values_matrix.csv", row.names = FALSE)
```

④ 合并文件
```r
library(dplyr)

median_data <- read.csv("./bp_hc_tsrna/median.csv")
de_values_data <- read.csv("./bp_hc_tsrna/de_values_matrix.csv")
# 基于 miRNA 列合并数据
merged_data <- inner_join(median_data, de_values_data, by = "tsRNA")

# 查看合并后的数据
print(merged_data)
write.csv(merged_data, file = "./bp_hc_tsrna/result.csv", row.names = FALSE)

```

```r
library(dplyr)
setwd("D:/adult_dep/DE2")

merged_data <- read.csv("./bp_hc_tsrna/result.csv")

diff_gene <- subset(merged_data, Q_value < 0.05 & abs(log2fold_change) > 0.58)
write.csv(diff_gene, file="./bp_hc_tsrna/sig_q0.05_fc15.csv", quote = F)
```



⑤ 绘图


```r
# 加载包
library(ggplot2)
library(pheatmap)

merged_data <- read.csv("./bp_hc_tsrna/result.csv")
merged_data$lgpadj <- -log10(merged_data$Q_value)

# Define categories for points
merged_data$category <- ifelse(merged_data$Q_value < 0.05 & merged_data$log2fold_change > 0.58, "Up", 
                                ifelse(merged_data$Q_value < 0.05 & merged_data$log2fold_change < -0.58, "Down", "Not Sig"))

# Calculate counts for each category
category_counts <- table(merged_data$category)

# Create a volcano plot
volcano_plot <- ggplot(merged_data, aes(x = log2fold_change, y = lgpadj)) +
  geom_point(aes(color = category), size = 2) +
  scale_color_manual(values = c("Up" = "#D32F2F", "Down" = "#1976D2", "Not Sig" = "grey70")) +
  labs(x = "Log2 Fold Change", y = "-Log10(padj)", title = "Volcano Plot") +
  theme_minimal(base_size = 14) +
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        axis.line = element_line(color = "black"),  # Show axis lines
        axis.ticks = element_line(color = "black"),  # Show axis ticks
        axis.ticks.length = unit(0.25, "cm")) +  # Set tick length
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", size = 0.5) +
  geom_vline(xintercept = c(-0.58, 0.58), linetype = "dashed", color = "black", size = 0.5) +
  
  # Add category count labels
  geom_text(aes(x = 2, y = 30, label = paste("Up:", category_counts["Up"])), color = "#D32F2F", size = 5, hjust = 0) +
  geom_text(aes(x = 2, y = 28, label = paste("Down:", category_counts["Down"])), color = "#1976D2", size = 5, hjust = 0) +
  geom_text(aes(x = 2, y = 26, label = paste("Not Sig:", category_counts["Not Sig"])), color = "grey70", size = 5, hjust = 0)

ggsave("./bp_hc_tsrna/bp_hc_tsrna_volcano_plot.pdf", plot = volcano_plot, device = "pdf", dpi = 600, width = 8, height = 6)
```




* 图中标注

```r
# 图中标注名称
annotations <- data.frame(
  tsRNA = c("ACTTGACCGCTCTGACCA", "ATCCTGCCGACTACGCCA", "GCATTGGTGGTTCAGTGGTAGAATTCT","GCCCGGATAGCTCAGTCGGTAGAG","GCGGGAGACCGGGGTTCGATTCCCCGACGGGGAGC","GTTTCCGTAGTGTAGTGGTTATC"),
  new_label = c("tsRNA-#9", "tsRNA-#17", "tsRNA-#36","tsRNA-#45","tsRNA-#51","tsRNA-#56")
)

# 提取需要标注的 miRNA 数据
annotations_data <- merged_data[merged_data$tsRNA %in% annotations$tsRNA, ]

# 创建火山图
volcano_plot <- ggplot(merged_data, aes(x = log2fold_change, y = lgpadj)) +
  geom_point(aes(color = ifelse(Q_value < 0.05 & log2fold_change > 0.58, "Up", 
                                ifelse(Q_value < 0.05 & log2fold_change < -0.58, "Down", "Not Sig"))), size = 2) +
  scale_color_manual(values = c("Up" = "#D32F2F", "Down" = "#1976D2", "Not Sig" = "grey70")) +
  labs(x = "Log2 Fold Change", y = "-Log10(padj)", title = "Volcano Plot") +
  theme_minimal(base_size = 14) +
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", size = 0.5) +
  geom_vline(xintercept = c(-0.58, 0.58), linetype = "dashed", color = "black", size = 0.5) +
  geom_point(data = annotations_data, aes(color = "Highlighted"), size = 3, shape = 1, stroke = 1.5, fill = NA) +  # 添加圈圈
  geom_text(data = merge(annotations_data, annotations, by = "tsRNA"), aes(label = new_label), 
            vjust = -1, hjust = 0.5, size = 3, color = "black")  # 使用新的标签


ggsave("./bp_hc_tsrna/label_bp_hc_tsrna_volcano_plot.pdf", plot = volcano_plot, device = "pdf", dpi = 600, width = 8, height = 6)

```
# 对mdd_uniq的mirna做富集分析

```r
library(grid)
library(futile.logger)
library(multiMiR)
library(tidyverse)

#设置路径
setwd("D:/adult_dep/DE2/intersect")
diff_gene <- read.csv("./mdd_uniq_mirna.csv", header = TRUE)


#up
miRNA_names_up <- diff_gene$miRNA
miRNA_suffixes_up <- sapply(strsplit(miRNA_names_up, "_"), function(x) tail(x, 1))

dg.miRNA.up <- get_multimir(org = "hsa", mirna = miRNA_suffixes_up, table = "mirecords", summary = TRUE)

table(dg.miRNA.up@data$type)  #178
up_result <- dg.miRNA.up@data
up_sum <- dg.miRNA.up@summary
write.csv(up_sum,file="./mdd_uniq_mirna_summary.csv")
write.csv(up_result,file="./mdd_uniq_mirna_result.csv")
```