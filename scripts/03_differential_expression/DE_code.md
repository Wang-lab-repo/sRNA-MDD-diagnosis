# Differential Expression Analysis Pipeline 

Take rsRNA in MDD vs. HC for example. 
  
## 1. Expression Quantification

```r
library(ggplot2)
setwd("./")
df <- read.csv("./input/rpm_50%_log_rsrna.csv", header = TRUE, row.names = 1)
df <- df[, df[1, ] %in% c(0, 1)]
df <- df[-1, ]
result_df <- data.frame()

for (sample_col in 2:ncol(df)) {
  sample_name <- colnames(df)[sample_col]
  non_zero_mirna_count <- sum(df[, sample_col] != 0)
  result_df <- rbind(result_df, data.frame(sample = sample_name, non_zero_mirna_count = non_zero_mirna_count))
}

write.csv(result_df, file = "./MDD_hc_rsrna/non_zero_rsrna_count.csv")

ggplot(result_df, aes(x = non_zero_mirna_count)) +
  geom_histogram(binwidth = 1, fill = "blue", color = "black", alpha = 0.7) +
  labs(title = "Frequency Distribution of Non-Zero rsRNA Counts",
       x = "Non-Zero rsRNA Count",
       y = "Frequency")
```

## 2. Sample Group Annotation
```r
df <- read.csv("./input/rpm_50%_log_rsrna.csv", header = TRUE, row.names = 1)
df <- df[, df[1, ] %in% c(0, 1)]
colnames(df) <- ifelse(df[1, ] == 0, paste0("control_", colnames(df)), paste0("treatment_", colnames(df)))
df <- df[-1, ]
write.csv(df, file = "./input_rename/rsrna.csv")
```

## 3. Median Calculation and Fold Change

```bash
python median_mirna.py ./input_rename/rsrna.csv ./MDD_hc_rsrna/median.csv
```

## 4. Two-sided Mann–Whitney U test
 
```r
data <- read.csv("./input_rename/rsrna.csv") 
data[data == 0] <- NA
write.csv(data, "./input_rename/rsrna_NA.csv", row.names = FALSE)  

data <- read.csv("./input_rename/rsrna_NA.csv")
library(dplyr)

control_columns <- grep("^control_", names(data), value = TRUE)
treatment_columns <- grep("^treatment_", names(data), value = TRUE)

p_values_matrix <- matrix(NA, nrow = length(data$rsRNA), ncol = 3)
colnames(p_values_matrix) <- c("rsRNA", "p_value", "Q_value")

for (i in seq_along(data$rsRNA)) {
  rsRNA <- data$rsRNA[i]
  control_values <- unlist(data[data$rsRNA == rsRNA, control_columns])
  treatment_values <- unlist(data[data$rsRNA == rsRNA, treatment_columns])
  result <- try(wilcox.test(treatment_values, control_values))
  if (!inherits(result, "try-error")) {
    p_values_matrix[i, 1] <- rsRNA
    p_values_matrix[i, 2] <- result$p.value
  }
}

p_values_matrix[, 3] <- p.adjust(p_values_matrix[, 2], method = "BH")
write.csv(p_values_matrix, file = "./MDD_hc_rsrna/de_values_matrix.csv", row.names = FALSE)
library(dplyr)
median_data <- read.csv("./MDD_hc_rsrna/median.csv")
de_values_data <- read.csv("./MDD_hc_rsrna/de_values_matrix.csv")
merged_data <- inner_join(median_data, de_values_data, by = "rsRNA")
write.csv(merged_data, file = "./MDD_hc_rsrna/result.csv", row.names = FALSE)

diff_gene <- subset(merged_data, Q_value < 0.05 & abs(log2fold_change) > 1)
write.csv(diff_gene, file="./MDD_hc_rsrna/sig_q0.05_fc2.csv", quote = F)

diff_gene <- subset(merged_data, Q_value < 0.05 & abs(log2fold_change) > 0.58)
write.csv(diff_gene, file="./MDD_hc_rsrna/sig_q0.05_fc15.csv", quote = F)
```


## 5. Volcano Plot


```r
library(ggplot2)
merged_data <- read.csv("./MDD_hc_rsrna/result.csv")
merged_data$lgpadj <- -log10(merged_data$Q_value)

ggplot(merged_data, aes(x = log2fold_change, y = lgpadj)) +
  geom_point(aes(color = ifelse(Q_value < 0.05 & log2fold_change > 0.58, "Up", 
                                ifelse(Q_value < 0.05 & log2fold_change < -0.58, "Down", "Not Sig"))), size = 2) +
  scale_color_manual(values = c("Up" = "#D32F2F", "Down" = "#1976D2", "Not Sig" = "grey70")) +
  labs(x = "Log2 Fold Change", y = "-Log10(padj)", title = "Volcano Plot") +
  theme_minimal() +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(-0.58, 0.58), linetype = "dashed", color = "black")

ggsave("./MDD_hc_rsrna/MDD_hc_rsrna_volcano_plot.pdf", dpi = 600, width = 8, height = 6)
```