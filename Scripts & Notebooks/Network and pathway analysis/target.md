# 0. miRNA靶点实验验证  
```r
BiocManager::install("grid", force = TRUE)
BiocManager::install("futile.logger")
BiocManager::install("multiMiR")
BiocManager::install("tidyverse")
#加载包
library(grid)
library(multiMiR)
library(tidyverse)

#设置路径
setwd("D:/adult_dep/DE2/mdd_common/所有预测结果")
cor.miRNA <- read.csv("./mirna.csv", header = TRUE)

miRNA_names <- cor.miRNA$miRNA
miRNA_suffixes <- sapply(strsplit(miRNA_names, "_"), function(x) tail(x, 1))

dg.miRNA <- get_multimir(org = "hsa", mirna = miRNA_suffixes, table = "mirecords", summary = TRUE)

table(dg.miRNA@data$type)  #1127
result <- dg.miRNA@data
sum <- dg.miRNA@summary
write.csv(sum,file="./only_mirna_实验/miRNA_summary.csv")
write.csv(result,file="./only_mirna_实验/miRNA_result.csv")
```


## 富集分析在差异分析的目录下完成，重新绘图富集分析
```r
library(ggplot2)
setwd("D:/adult_dep/HAMD_corr")
# 从 CSV 文件读取数据
data <- read.csv("./mirna/go_top10.csv", stringsAsFactors = FALSE)
# 计算 -log10(P value)
data$neg_log10_pvalue <- -log10(data$P.value)

# data <- data %>%
#   arrange(desc(neg_log10_pvalue))


# 定义渐变颜色
start_color <- rgb(173/255, 216/255, 230/255)  
end_color <- rgb(0/255, 51/255, 102/255)  

# 绘制条形图
ggplot(data, aes(x = E.ratio, y = reorder(GO.term, E.ratio), fill = neg_log10_pvalue)) +
  geom_bar(stat = "identity", width = 0.8) +  # 调整条形的宽度
  scale_fill_gradient(low = start_color, high = end_color, name = "-log10(P value)") +
  labs(x = "E-ratio", y = "Pathway", title = "Enrichment Analysis Bar Plot") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 10),
    panel.grid.major = element_blank(),  # 去除主要网格线
    panel.grid.minor = element_blank(),  # 去除次要网格线
    panel.background = element_rect(fill = "white")  # 背景颜色设置为白色
  )

ggsave("./mirna/cor_go_enrichment_barplot_top10.svg", plot = last_plot(), device = "svg", dpi = 600)
```



# 1. sRNA靶点预测
## miRanda
* 整理格式

```py
import re

def clean_fasta(input_file, output_file):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        header = None
        sequence = ''
        for line in infile:
            line = line.strip()
            if line.startswith('>'):
                if header and sequence:
                    outfile.write(header + '\n' + sequence + '\n')
                header = line
                sequence = ''
            else:
                sequence += re.sub(r'[^ATCGUN]', '', line)
        if header and sequence:  # Write the last entry
            outfile.write(header + '\n' + sequence + '\n')

# 清理 FASTA 文件
clean_fasta('mRNA_UTR.fasta', 'mRNA_UTR_clean.fasta')

python clean_utr.py
```

```bash
conda activate python39
# conda install bioconda::miranda
cd /mnt/d/adult_dep/DE2/mdd_common/所有预测结果/input
cat mirna_seq.fa rsrna.fa tsrna.fa pirna_seq.fa > seq.fa

ref_dir=/mnt/d/adult_dep/final_penal_target
#预测靶点
miranda seq.fa ${ref_dir}/mRNA_UTR_clean.fasta -out out.txt 
# -sc 150 -en -20

#输出结果
grep '>>' out.txt > miranda_result.txt
```

```py
vim run.py
####################################

#整理格式
import os
import csv

def TidyMirandaResult(path, inputfile, outfile):
    infpath = os.path.join(path, inputfile)
    outfpath = os.path.join(path, outfile)
    
    with open(infpath, 'r') as f:
        lines = f.readlines()
        
    with open(outfpath, 'w', newline='') as csvfile:
        csv_writer = csv.writer(csvfile)
        csv_writer.writerow(['sRNAs', 'Gene', 'Score', 'Energy'])  # CSV文件的列头
        
        for line in lines:
            line = line.strip()
            if line.startswith('>>'):
                parts = line.split('\t')
                
                # 处理并提取所需数据
                srna_name = parts[0][2:]  # 去除'>>'前缀
                gene_name = parts[1].split('|')[0]  # 提取基因的第一个标识符
                score = parts[2]  # 假设是评分字段
                energy = parts[3]  # 假设是能量字段
                
                # 写入CSV文件
                csv_writer.writerow([srna_name, gene_name, score, energy])

def main():
    path = r'/mnt/d/adult_dep/final_penal_target'
    input_file = 'miranda_result.txt'
    output_file = 'miranda_TidyResult.csv'
    TidyMirandaResult(path, input_file, output_file)

if __name__ == "__main__":
    main()
####################################


python run.py
```
## Targetscan
* 用于miRNA的预测

## RNAhybrid
```bash
-b <number of hits per target>: 指定每个目标序列的最大命中数。默认情况下，它可能会返回多个结合位置。
-c: 输出紧凑格式。会以更简洁的格式显示预测结果。
-d <xi>,<theta>: 用于指定极值分布的参数 xi 和 theta，用于计算p值。如果没有指定，工具将根据查询序列的最大配对能量来估计这些值。
-f <helix constraint>: 指定约束条件，定义一个特定的螺旋结构（例如：-f 2,7 强制要求查询序列在位置2到7之间有一个螺旋结构）。
-h: 显示帮助信息，列出所有可用的选项。
-m <max targetlength>: 设置目标序列的最大长度。
-n <max query length>: 设置查询序列的最大长度。
-u <max internal loop size (per side)>: 设置最大内环的大小（每一侧的大小）。
-v <max bulge loop size>: 设置最大突起环（bulge loop）的大小。
-e <energy cut-off>: 设置能量截断值，用于过滤不太可能发生的配对。
-p <p-value cut-off>: 设置p值截断阈值，用于过滤结果。
-s (3utr_fly|3utr_worm|3utr_human): 指定数据集的名称（例如：3utr_fly，3utr_worm，3utr_human），用于估计p值的分布参数。
-g (ps|png|jpg|all): 指定图形输出格式。ps 是PostScript格式，png 和 jpg 是图像格式，all 会输出所有格式。
-t <target file>: 指定一个包含目标序列的文件（FASTA格式）。
-q <query file>: 指定一个包含查询序列的文件（FASTA格式）。
```

```bash
# conda install bioconda::RNAhybrid 
cd /mnt/d/adult_dep/DE2/mdd_common/所有预测结果/input
ref_dir=/mnt/d/adult_dep/final_penal_target

RNAhybrid -g jpg -b 1 -e -20 -f 8,12 -u 1 -v 1 -s 3utr_human -t ${ref_dir}/mRNA_UTR_clean.fasta -q seq.fa  > hybird_result.txt
```

```py
vim RNAhybrid.py
####################################
# 整理格式
import re
import pandas as pd

def parse_rnahybrid_results(result_file):
    results = []
    current_srna = None
    current_target = None
    current_pvalue = None
    
    with open(result_file, 'r') as file:
        for line in file:
            line = line.strip()
            
            if line.startswith("miRNA :"):
                # 提取miRNA信息
                current_srna = re.search(r'^miRNA : (\S+)', line).group(1)
            
            elif line.startswith("target:"):
                # 提取target信息，保留ENSG编号
                target_line = re.search(r'^target: (\S+)', line).group(1)
                current_target = target_line.split('|')[0]  # 只保留第一个ENSG编号
            
            elif line.startswith("p-value:"):
                # 提取p-value信息
                current_pvalue = float(re.search(r'p-value: (\S+)', line).group(1))
                
                # 如果p-value小于0.05，保存当前的miRNA和target
                if current_pvalue and current_pvalue < 0.05:
                    results.append([current_srna, current_target, current_pvalue])
    
    return results

# 解析hybird_result.txt文件中的靶基因信息
rsRNA_targets = parse_rnahybrid_results('hybird_result.txt')

# 保存去重后的靶基因列表到文件
pd.DataFrame(rsRNA_targets, columns=['sRNAs', 'Gene', 'p-value']).drop_duplicates().to_csv("hybrid_result.csv", index=False)
####################################
python RNAhybrid.py
```
## 取交集
```py
import pandas as pd

# 读取 CSV 文件
hybrid_result_df = pd.read_csv('hybrid_result.csv')
miranda_result_df = pd.read_csv('miranda_TidyResult.csv')

# 筛选 Energy 小于 -20 的行
miranda_result_df = miranda_result_df[miranda_result_df['Energy'] < -20]

# 按照'sRNAs'列分组，取交集
result = []

# 获取所有独特的 sRNAs
sRNAs_list = hybrid_result_df['sRNAs'].unique()

# 遍历所有sRNAs，找交集
for srna in sRNAs_list:
    # 获取该sRNAs在hybrid_result_df中的Gene
    hybrid_genes = set(hybrid_result_df[hybrid_result_df['sRNAs'] == srna]['Gene'])
    
    # 获取该sRNAs在miranda_result_df中的Gene
    miranda_genes = set(miranda_result_df[miranda_result_df['sRNAs'] == srna]['Gene'])
    
    # 计算交集
    intersect_genes = hybrid_genes.intersection(miranda_genes)
    
    # 如果交集不为空，保存结果
    if intersect_genes:
        result.append([srna, list(intersect_genes)])

# 保存交集结果到 CSV 文件
output_df = pd.DataFrame(result, columns=['sRNAs', 'Gene_Intersection'])
output_df.to_csv('gene_intersection.csv', index=False)


####################################### 另一种保存格式  #######################################
import pandas as pd

# 读取 CSV 文件
hybrid_result_df = pd.read_csv('hybrid_result.csv')
miranda_result_df = pd.read_csv('miranda_TidyResult.csv')

# 筛选 Energy 小于 -20 的行
miranda_result_df = miranda_result_df[miranda_result_df['Energy'] < -20]

# 获取所有独特的 sRNAs
sRNAs_list = hybrid_result_df['sRNAs'].unique()

# 用来保存结果的列表
result = []

# 遍历所有sRNAs，找交集
for srna in sRNAs_list:
    # 获取该sRNAs在hybrid_result_df中的Gene
    hybrid_genes = set(hybrid_result_df[hybrid_result_df['sRNAs'] == srna]['Gene'])
    
    # 获取该sRNAs在miranda_result_df中的Gene
    miranda_genes = set(miranda_result_df[miranda_result_df['sRNAs'] == srna]['Gene'])
    
    # 计算交集
    intersect_genes = hybrid_genes.intersection(miranda_genes)
    
    # 如果交集不为空，保存每个基因
    for gene in intersect_genes:
        result.append([srna, gene])

# 将结果转换为 DataFrame
result_df = pd.DataFrame(result, columns=['sRNA', 'Gene'])

result_df.to_csv('sRNA_gene_intersection.csv', index=False)

# 显示结果（或调试）
print(result_df)
```
## ENSEMBL转换为symbol名称

```r
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(clusterProfiler)

setwd("D:/adult_dep/DE2/mdd_common/所有预测结果")
srna <- read.csv("./sRNA_gene_intersection.csv")

ensembl_id_transform <- function(ENSEMBL_ID) {
    a = bitr(ENSEMBL_ID, fromType = "ENSEMBL", toType = c("SYMBOL", "ENTREZID"), OrgDb = "org.Hs.eg.db")
    return(a)
  }

transformed <- ensembl_id_transform(srna$Gene)

# 确保合并后的数据框与原始数据框正确对齐
srna <- merge(srna, transformed, by.x = "Gene", by.y = "ENSEMBL", all.x = TRUE)

# 查看结果
head(srna)
write.csv(srna, file = "sRNA_gene_intersection_transname.csv", quote = FALSE)
# 删除NA
```







## 筛选想要的前10通路，重新绘图富集分析
```r
library(ggplot2)
setwd("D:/adult_dep/DE2/mdd_common/所有预测结果/only_mirna_实验")
# 从 CSV 文件读取数据
data <- read.csv("./go.csv", stringsAsFactors = FALSE)

# 计算 -log10(P value)
data$neg_log10_pvalue <- -log10(data$P.value)

# 定义渐变颜色
start_color <- rgb(173/255, 216/255, 230/255)  # 深蓝色
end_color <- rgb(0/255, 51/255, 102/255)    # 橘黄色

# 绘制条形图
ggplot(data, aes(x = E.ratio, y = reorder(GO.term, E.ratio), fill = neg_log10_pvalue)) +
  geom_bar(stat = "identity", width = 0.8) +  # 调整条形的宽度
  scale_fill_gradient(low = start_color, high = end_color, name = "-log10(P value)") +
  labs(x = "Gene ratio", y = "Pathway", title = "Enrichment Analysis Bar Plot") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 10),
    panel.grid.major = element_blank(),  # 去除主要网格线
    panel.grid.minor = element_blank(),  # 去除次要网格线
    panel.background = element_rect(fill = "white")  # 背景颜色设置为白色
  )

ggsave("./only_mirna_实验_enrichment_barplot.svg", plot = last_plot(), device = "svg", dpi = 600)
```

```r
library(ggplot2)
setwd("E:/adolescent_dep/target/all")
# 从 CSV 文件读取数据
data <- read.csv("./mdd_bp.csv", stringsAsFactors = FALSE)

# 计算 -log10(P value)
data$neg_log10_pvalue <- -log10(data$P.value)

# 定义渐变颜色
start_color <- rgb(173/255, 216/255, 230/255)  # 深蓝色
end_color <- rgb(0/255, 51/255, 102/255)    # 橘黄色

# 绘制条形图
ggplot(data, aes(x = E.ratio, y = reorder(GO.term, E.ratio), fill = neg_log10_pvalue)) +
  geom_bar(stat = "identity", width = 0.8) +  # 调整条形的宽度
  scale_fill_gradient(low = start_color, high = end_color, name = "-log10(P value)") +
  labs(x = "E-ratio", y = "Pathway", title = "Enrichment Analysis Bar Plot") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 10),
    panel.grid.major = element_blank(),  # 去除主要网格线
    panel.grid.minor = element_blank(),  # 去除次要网格线
    panel.background = element_rect(fill = "white"))


ggsave("./bp_enrichment_barplot.svg", 
       plot = last_plot(), 
       device = "svg", 
       dpi = 600, 
       width = 12,    # 宽度调整为 12 英寸
       height = 8,    # 高度调整为 8 英寸
       units = "in")  # 单位设置为英寸

# 以p排序
# Ensure your data is arranged by neg_log10_pvalue in descending order
data <- data %>%
  arrange(desc(neg_log10_pvalue))

# Plot with ggplot2
ggplot(data, aes(x = E.ratio, y = reorder(GO.term, neg_log10_pvalue), fill = neg_log10_pvalue)) +
  geom_bar(stat = "identity", width = 0.8) +  # Adjust the bar width
  scale_fill_gradient(low = start_color, high = end_color, name = "-log10(P value)") +
  labs(x = "E-ratio", y = "Pathway", title = "Enrichment Analysis Bar Plot") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 10),
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.background = element_rect(fill = "white")  # Set background color to white
  )


ggsave("./MDD_hc_mirna/down_final/go_enrichment_barplot_top10.svg", plot = last_plot(), device = "svg", dpi = 600)
```


## KEGG分析
```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "enrichplot"))  # 人类数据库示例

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

setwd("D:/adult_dep/DE2/mdd_common/所有预测结果/")
gene_list <- read.csv("sRNA_gene_intersection_transname.csv", header=TRUE)$ENTREZID  # 替换为实际路径
gene_entrez <- unique(gene_list)

kegg_result <- enrichKEGG(gene = gene_entrez,
                          organism = "hsa",  # 人类：hsa，小鼠：mmu
                          keyType = "kegg",
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.05)

# 查看结果（按p值排序）
head(kegg_result)
write.csv(kegg_result, "./target/kegg_results.csv", row.names = FALSE)


barplot(kegg_result, showCategory=20)                # 条形图
dotplot(kegg_result, showCategory=10, color="pvalue", size="Count")                # 点图

p_dot <- dotplot(kegg_result, showCategory=10, color="p.adjust", size="GeneRatio") +
  scale_size_continuous(range=c(3,8)) +
  theme(axis.text.y=element_text(size=10))
ggsave("kegg_dotplot.pdf", p_dot, width=8, height=6)


emapplot(pairwise_termsim(kegg_result))              # 网络图
cnetplot(kegg_result, categorySize="pvalue", foldChange=gene_entrez)  # 基因-通路网络图



# 修改后代码
library(clusterProfiler)
library(ggplot2)

# 假设 kegg_result 是已有的 enrichKEGG 结果
# 计算 Fold Enrichment
# 计算 Fold Enrichment
kegg_result@result$GeneRatio_num <- sapply(strsplit(kegg_result@result$GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))
kegg_result@result$BgRatio_num <- sapply(strsplit(kegg_result@result$BgRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))
kegg_result@result$FoldEnrichment <- kegg_result@result$GeneRatio_num / kegg_result@result$BgRatio_num
# kegg_result@result$pvalue <- -log10(kegg_result@result$pvalue)  # 新增-log10转换

# 查看计算结果
head(kegg_result@result[, c("Description", "GeneRatio", "BgRatio", "FoldEnrichment")])

# 备份原始 GeneRatio
original_gene_ratio <- kegg_result@result$GeneRatio

# 替换 GeneRatio 为 FoldEnrichment（使 dotplot 使用它）
kegg_result@result$GeneRatio <- kegg_result@result$FoldEnrichment


# 绘制点图
p <- dotplot(kegg_result, 
             showCategory = 10, 
             color = "pvalue", 
             size = "Count",
             x = "FoldEnrichment") +  # 使用修改后的 GeneRatio（实际是 FoldEnrichment）
  labs(x = "Enrichment ratio",     # 修改 x 轴标签
       color = "-log10(P value)",         # 修改图例标签
       size = "Gene count") +     # 修改图例标签
  theme_minimal(base_size = 12)   # 调整字体大小和主题

# 保存图片
ggsave("kegg_dotplot_fold_enrichment.pdf", p, width = 8, height = 6)
```
















# 3.  STRINGdb 聚类分析
STRING是一个已知和预测的蛋白质-蛋白质相互作用的数据库。 相互作用包括直接(物理)和间接(功能)联系;它们源于计算预测、生物之间的知识转移，以及其他(主要)数据库聚合的交互作用。 STRING中的相互作用有五个主要来源：基因组预测、高通量实验、（保守的）共表达实验、自动化文本挖掘、数据库相关知识。    
string_db$get_clusters() 返回的是基因的 社区划分信息，指示每个基因属于哪个社区。  
string_db$plot_network() 绘制的是基因间的 相互作用网络图，展示基因之间的相互作用，但没有明确显示社区信息。  

## MDD10_score_threshold=700
```r
BiocManager::install("STRINGdb")

library(igraph)
library(STRINGdb)
library(tidyverse)

setwd("D:/adult_dep/penal_target")
degs <- read.table("./10/mdd10.txt", header=FALSE, sep="\t")

head(degs)

# 选择载入的STRINGdb数据（数据库版本、物种、蛋白互作得分）#人9606，小鼠10090 #蛋白互作的得分 默认400, 低150，高700，极高900
string_db <- STRINGdb$new(version="12.0", species=9606, score_threshold=700, input_directory="")

# 选择基因ID类型，如果是Ensembl ID或者基因符号都可以
gene_ids <- degs$V1

# 将基因ID映射到STRING数据库
dat_map <- string_db$map(my_data_frame=degs, 
                         my_data_frame_id_col_names="V1", #使用gene symbol或ENTREZID都可
                         removeUnmappedRows = TRUE )
# dat_map <- string_db$map(my_data_frame=gene_ids, 
#                          my_data_frame_id_col_names="gene", #使用gene symbol或ENTREZID都可
#                          removeUnmappedRows = TRUE )
hits <- dat_map$STRING_id
length(hits) #448

# 使用plot_network即可绘制PPI图，还可以给PPI添加上下调信息（ green down-regulated gened and red for up-regulated genes）
## PPI
png("string_PPI_mdd10.png",units="in",width = 10,height = 10, res=400)
string_db$plot_network(hits)
dev.off()
## PPI_halo  #给PPI添加上下调信息
# filter by p-value and add a color column(i.e.green for down and red for up genes)
dat_map_color <- string_db$add_diff_exp_color(subset(dat_map, pvalue<0.05),
                                                  logFcColStr="log2FC" )
payload_id <- string_db$post_payload(dat_map_color$STRING_id,
                                     colors=dat_map_color$color)
png("string_PPI_halo.png",units="in",width = 10,height = 10, res=400)
string_db$plot_network(hits, payload_id=payload_id )
dev.off()




## iGraph clustering 互作网络分簇
# algorithm: fastgreedy(默认), walktrap, edge.betweenness
# 使用 Fast Greedy 方法进行聚类，但是结果太分散，因此尝试其他算法
clustersList <- string_db$get_clusters(string_ids = hits, algorithm  = "fastgreedy" ) 
png("string_PPI_iGraph_cluster_mdd10.png",units="in",width = 15,height = 10,res=400)
par(mfrow=c(2,3))
for(i in 1:6){
 string_db$plot_network(clustersList[[i]])
}
dev.off()




# 使用 walktrap 方法进行聚类——没啥区别
network_data <- string_db$get_interactions(dat_map$STRING_id)
head(network_data)

library(igraph)
igraph_network <- graph_from_data_frame(network_data)
plot(igraph_network)

walktrap_clusters <- cluster_walktrap(igraph_network)
membership(walktrap_clusters)
community_genes <- membership(walktrap_clusters)
for (cluster_id in unique(community_genes)) {
  cluster_genes_list <- names(community_genes[community_genes == cluster_id])
  print(paste("Cluster", cluster_id, "genes:", paste(cluster_genes_list, collapse=", ")))
}
png("string_PPI_iGraph_cluster_mdd10——2.png",units="in",width = 15,height = 10,res=400)
par(mfrow=c(2,3))
for(i in 1:6){
 string_db$plot_network(walktrap_clusters[[i]])
}
dev.off()


############################## 获取蛋白互作信息用于后续可视化 ###############
dat_link <- string_db$get_interactions(hits)  # dat_link=network_data
# 转换stringID为 gene symbol
dat_link$from <- dat_map[match(dat_link$from,dat_map$STRING_id),'V1']
dat_link$to <- dat_map[match(dat_link$to,dat_map$STRING_id),'V1']  
colnames(dat_link) <- c('node1','node2','combined_score')
# 去除重复
dat_link <- dat_link %>% distinct(node1, node2, .keep_all = T)

write.csv(dat_link,'string_link_mdd10.csv',row.names = F,quote = F)

```

## MDD10_score_threshold=400
```r
BiocManager::install("STRINGdb")

library(igraph)
library(STRINGdb)
library(tidyverse)

setwd("D:/adult_dep/penal_target")
degs <- read.table("./10/mdd10.txt", header=FALSE, sep="\t")

string_db <- STRINGdb$new(version="12.0", species=9606, score_threshold=400, input_directory="")

gene_ids <- degs$V1

# 将基因ID映射到STRING数据库
dat_map <- string_db$map(my_data_frame=degs, 
                         my_data_frame_id_col_names="V1", #使用gene symbol或ENTREZID都可
                         removeUnmappedRows = TRUE )

hits <- dat_map$STRING_id
length(hits) #448

png("string_PPI_mdd10.png",units="in",width = 10,height = 10, res=400)
string_db$plot_network(hits)
dev.off()

# dat_map_color <- string_db$add_diff_exp_color(subset(dat_map, pvalue<0.05),
#                                                   logFcColStr="log2FC" )
# payload_id <- string_db$post_payload(dat_map_color$STRING_id,
#                                      colors=dat_map_color$color)
# png("string_PPI_halo.png",units="in",width = 10,height = 10, res=400)
# string_db$plot_network(hits, payload_id=payload_id )
# dev.off()


## iGraph clustering 互作网络分簇
clustersList <- string_db$get_clusters(string_ids = hits, algorithm  = "fastgreedy" ) 
png("string_PPI_iGraph_cluster_mdd10.png",units="in",width = 15,height = 10,res=400)
par(mfrow=c(2,3))
for(i in 1:6){
 string_db$plot_network(clustersList[[i]])
}
dev.off()

############################## 获取蛋白互作信息用于后续可视化 ###############
dat_link <- string_db$get_interactions(hits)  # dat_link=network_data
# 转换stringID为 gene symbol
dat_link$from <- dat_map[match(dat_link$from,dat_map$STRING_id),'V1']
dat_link$to <- dat_map[match(dat_link$to,dat_map$STRING_id),'V1']  
colnames(dat_link) <- c('node1','node2','combined_score')
# 去除重复
dat_link <- dat_link %>% distinct(node1, node2, .keep_all = T)

write.csv(dat_link,'string_link_mdd10.csv',row.names = F,quote = F)

```

## Other10_score_threshold=400

```r
BiocManager::install("STRINGdb")

library(igraph)
library(STRINGdb)
library(tidyverse)

setwd("D:/adult_dep/penal_target")
degs <- read.table("./10/other10.txt", header=FALSE, sep="\t")

string_db <- STRINGdb$new(version="12.0", species=9606, score_threshold=400, input_directory="")

gene_ids <- degs$V1

# 将基因ID映射到STRING数据库
dat_map <- string_db$map(my_data_frame=degs, 
                         my_data_frame_id_col_names="V1", #使用gene symbol或ENTREZID都可
                         removeUnmappedRows = TRUE )

hits <- dat_map$STRING_id
length(hits) #19

png("string_PPI_other10.png",units="in",width = 10,height = 10, res=400)
string_db$plot_network(hits)
dev.off()

# dat_map_color <- string_db$add_diff_exp_color(subset(dat_map, pvalue<0.05),
#                                                   logFcColStr="log2FC" )
# payload_id <- string_db$post_payload(dat_map_color$STRING_id,
#                                      colors=dat_map_color$color)
# png("string_PPI_halo.png",units="in",width = 10,height = 10, res=400)
# string_db$plot_network(hits, payload_id=payload_id )
# dev.off()


## iGraph clustering 互作网络分簇
clustersList <- string_db$get_clusters(string_ids = hits, algorithm  = "fastgreedy" ) 
png("string_PPI_iGraph_cluster.png",units="in",width = 15,height = 10,res=400)
par(mfrow=c(2,3))
for(i in 1:6){
 string_db$plot_network(clustersList[[i]])
}
dev.off()

############################## 获取蛋白互作信息用于后续可视化 ###############
dat_link <- string_db$get_interactions(hits)  # dat_link=network_data
# 转换stringID为 gene symbol
dat_link$from <- dat_map[match(dat_link$from,dat_map$STRING_id),'V1']
dat_link$to <- dat_map[match(dat_link$to,dat_map$STRING_id),'V1']  
colnames(dat_link) <- c('node1','node2','combined_score')
# 去除重复
dat_link <- dat_link %>% distinct(node1, node2, .keep_all = T)

write.csv(dat_link,'string_link.csv',row.names = F,quote = F)
```

## all_mirna_score_threshold=400

```r
BiocManager::install("STRINGdb")

library(igraph)
library(STRINGdb)
library(tidyverse)

setwd("D:/adult_dep/penal_target")
degs <- read.table("./all/all.txt", header=FALSE, sep="\t")

string_db <- STRINGdb$new(version="11.5", species=9606, score_threshold=400, input_directory="")

gene_ids <- degs$V1

# 将基因ID映射到STRING数据库
dat_map <- string_db$map(my_data_frame=degs, 
                         my_data_frame_id_col_names="V1", #使用gene symbol或ENTREZID都可
                         removeUnmappedRows = TRUE )

hits <- dat_map$STRING_id #62%gene 没有对应的id
length(hits) #13

png("./all/string_PPI_all10.png",units="in",width = 10,height = 10, res=400)
string_db$plot_network(hits)
dev.off()

# dat_map_color <- string_db$add_diff_exp_color(subset(dat_map, pvalue<0.05),
#                                                   logFcColStr="log2FC" )
# payload_id <- string_db$post_payload(dat_map_color$STRING_id,
#                                      colors=dat_map_color$color)
# png("string_PPI_halo.png",units="in",width = 10,height = 10, res=400)
# string_db$plot_network(hits, payload_id=payload_id )
# dev.off()


# ## iGraph clustering 互作网络分簇
# clustersList <- string_db$get_clusters(string_ids = hits, algorithm  = "fastgreedy" ) 
# png("string_PPI_iGraph_cluster.png",units="in",width = 15,height = 10,res=400)
# par(mfrow=c(2,3))
# for(i in 1:6){
#  string_db$plot_network(clustersList[[i]])
# }
# dev.off()

############################## 获取蛋白互作信息用于后续可视化 ###############
dat_link <- string_db$get_interactions(hits)  # dat_link=network_data
# 转换stringID为 gene symbol
dat_link$from <- dat_map[match(dat_link$from,dat_map$STRING_id),'V1']
dat_link$to <- dat_map[match(dat_link$to,dat_map$STRING_id),'V1']  
colnames(dat_link) <- c('node1','node2','combined_score')
# 去除重复
dat_link <- dat_link %>% distinct(node1, node2, .keep_all = T)

write.csv(dat_link,'./all/string_link.csv',row.names = F,quote = F)
```