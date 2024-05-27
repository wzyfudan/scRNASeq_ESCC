rm(list=ls())
setwd("E:/Lab data/Sunlab/bioinformatics/rawdata")
library(Seurat)
pbmc <- readRDS("Mono_cell.rds")
table(pbmc$new.cluster.ids)
#使用R从Seurat对象中提取标准差最大的top基因计算相关性，用pheatmap画图。
av <-AverageExpression(pbmc,
                       group.by = "new.cluster.ids",
                       assays = "RNA")
av=av[[1]]
head(av)

#选出标准差最大的1000个基因
cg=names(tail(sort(apply(av, 1, sd)),1000))
#查看这1000个基因在各细胞群中的表达矩阵
View(av[cg,])
#查看细胞群的表达相关性矩阵
View(cor(av[cg,],method = 'spearman'))
#pheatmap绘制热图
pheatmap::pheatmap(cor(av[cg,],method = 'spearman')) #默认是Pearson

#细胞群共现矩阵Co-occurrence-------------------------------------------------
library(dplyr)
ESOtest <- readRDS("Mono_cell.Rds")
ESOtest <- readRDS("CD8T_cell.Rds")
ESOtest <- readRDS("CD4T_cell.Rds")
ESOtest <- readRDS("B_cell_harmony.Rds")
#各类群细胞比例
meta = ESOtest@meta.data
NR = rbind(meta[which(meta$orig.ident == 'P1'),],meta[which(meta$orig.ident == 'P5'),],meta[which(meta$orig.ident == 'P7'),],meta[which(meta$orig.ident == 'P9'),],meta[which(meta$orig.ident == 'P10'),])
R = rbind(meta[which(meta$orig.ident == 'P2'),],meta[which(meta$orig.ident == 'P3'),],meta[which(meta$orig.ident == 'P4'),],meta[which(meta$orig.ident == 'P6'),],meta[which(meta$orig.ident == 'P8'),])
ids1 = table(NR[,c("new.cluster.ids","orig.ident")]) %>% as.data.frame() 
ids2 = table(R[,c("new.cluster.ids","orig.ident")]) %>% as.data.frame()
ids1$group = 'NR'
ids2$group = 'R'
df = rbind(ids1,ids2)

df1 = rbind(ids1,ids2) #换个interest细胞亚群重复上述比例计算
df <- rbind(df,df1) #如合并CD8T和B细胞的Freq
df$Freq[df$Freq == 0]<-NA
df = na.omit(df)

#计算prop
ESOtest1 <- readRDS("scRNA_harmony_UMAP.Rds")
meta1 = ESOtest1@meta.data
ids = table(meta1$orig.ident) %>% as.data.frame()
for(i in 1:dim(ids)[1]){
  prop.temp = df[which(df$orig.ident %in% ids[i,1]),'Freq'] / ids[i,2]
  df[which(df$orig.ident %in% ids[i,1]),'prop'] = prop.temp
}
#制作矩阵
library(tidyr)
library(dplyr)

prop_matrix <- df %>%
  # 将表格按照Patient_ID和Cell_Type分组，然后对Prop取平均值
  group_by(new.cluster.ids,orig.ident) %>%
  summarise(Avg_Prop = mean(prop)) %>%
  # 使用pivot_wider函数将数据从长格式转换为宽格式，形成prop矩阵
  pivot_wider(names_from = new.cluster.ids, values_from = Avg_Prop)

prop_matrix <- as.data.frame(prop_matrix)
prop_matrix <- prop_matrix[,-1]
prop_matrix[is.na(prop_matrix)] <- 0
prop_matrix <- as.matrix(prop_matrix)

spearman_corr <- cor(prop_matrix, method = "spearman")
subset_corr <- spearman_corr[10:20,c("TNFRSF13B Bmem","FCRL4 Bmem"), drop = FALSE]
#可视化
library(ggplot2)
library(RColorBrewer)
pheatmap::pheatmap(subset_corr, 
                   col = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),  # 使用调色板设置颜色
                   scale = "none",  # 不进行数据标准化
                   main = "Bmem Spearman Correlation",
                   fontsize_col = 13,   # Set column names font size
                   angle_col = 45,# 设置列名方向（角度）
                   cellwidth = 30,# 设置每个单元格的宽度
                   cellheight = 25,
                   cluster_rows = FALSE,          # 禁用行聚类
                   cluster_cols = FALSE )         # 禁用列聚类 
pdf("CD8-DC Rho.PDF", width = 6,height = 4)        
