library(ggplot2)
setwd("./")
data <- read.csv("./go.csv", stringsAsFactors = FALSE)
data$neg_log10_pvalue <- -log10(data$P.value)

start_color <- rgb(173/255, 216/255, 230/255)
end_color <- rgb(0/255, 51/255, 102/255)

ggplot(data, aes(x = E.ratio, y = reorder(GO.term, E.ratio), fill = neg_log10_pvalue)) +
  geom_bar(stat = "identity", width = 0.8) +
  scale_fill_gradient(low = start_color, high = end_color, name = "-log10(P value)") +
  labs(x = "Gene ratio", y = "Pathway", title = "Enrichment Analysis Bar Plot") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 10),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white")
  )

ggsave("enrichment_barplot.svg", plot = last_plot(), device = "svg", dpi = 600)