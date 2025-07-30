# Functional Enrichment Analysis Workflow

This document describes the workflow for functional enrichment analysis using R, including gene ID conversion, GO and KEGG enrichment, and visualization with bar plots and dot plots.  


## 1. Convert ENSEMBL IDs to Gene Symbols

```r
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(clusterProfiler)

setwd("./")
srna <- read.csv("./sRNA_gene_intersection.csv")

ensembl_id_transform <- function(ENSEMBL_ID) {
    bitr(ENSEMBL_ID, fromType = "ENSEMBL", toType = c("SYMBOL", "ENTREZID"), OrgDb = org.Hs.eg.db)
}

transformed <- ensembl_id_transform(srna$Gene)
srna <- merge(srna, transformed, by.x = "Gene", by.y = "ENSEMBL", all.x = TRUE)
write.csv(srna, file = "sRNA_gene_intersection_transname.csv", quote = FALSE)
```

## 2. GO Pathway Enrichment Analysis
Extract the Gene Symbols from 'sRNA_gene_intersection_transname.csv' and save them as a .txt file.  

```bash
for i in *.txt;do
  python enrichgo.py goa_human.gaf go.obo uniprot-proteome_UP000005640_reviewed_yes.fasta ${i} ${i%%.*}_go.txt
done
```



## 3. KEGG Pathway Enrichment Analysis
```r
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)


setwd("./")
gene_list <- read.csv("sRNA_gene_intersection_transname.csv", header=TRUE)$ENTREZID
gene_entrez <- unique(gene_list)

kegg_result <- enrichKEGG(gene = gene_entrez,
                          organism = "hsa",
                          keyType = "kegg",
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.05)

write.csv(kegg_result, "kegg_results.csv", row.names = FALSE)

barplot(kegg_result, showCategory=20)
dotplot(kegg_result, showCategory=10, color="pvalue", size="Count")
```
