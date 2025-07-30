library(ggplot2)
library(pheatmap)

setwd("./")
data <- read.csv("spearman_results_corrected.csv", header = TRUE, row.names = 1)
data$lgadjp <- -log10(data$Adjusted_p_value)

volcano_plot <- ggplot(data, aes(x = Correlation, y = lgadjp)) +
  geom_point(aes(color = ifelse(Adjusted_p_value < 0.05 & Correlation > 0.2, "Pos", 
                                ifelse(Adjusted_p_value < 0.05 & Correlation < -0.2, "Neg", "Not Sig"))), size = 2) +
  scale_color_manual(values = c("Pos" = "#D32F2F", "Neg" = "#1976D2", "Not Sig" = "grey70")) +
  labs(x = "Spearman correlation", y = "-lg(adj)", title = "Volcano Plot") +
  theme_minimal(base_size = 14) +
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(), 
        axis.line = element_line(color = "black"), 
        axis.ticks = element_line(color = "black"),  
        axis.ticks.length = unit(0.25, "cm")) +  
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", size = 0.5) +
  geom_vline(xintercept = c(-0.2, 0.2), linetype = "dashed", color = "black", size = 0.5)

ggsave("volcano_plot.pdf", plot = volcano_plot, device = "pdf", dpi = 600, width = 8, height = 6)