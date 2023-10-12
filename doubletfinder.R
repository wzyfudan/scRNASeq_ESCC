#https://cloud.tencent.com/developer/article/1825672
library(Seurat)
library(patchwork)
library(tidyverse)
devtools::install_github('chris-mcginnis-ucsf/DoubletFinder')
library(DoubletFinder)
library(KernSmooth)
library(ROCR)
setwd("E:/SHMC_1805/Lab data/Sunlab/bioinformatics/rawdata")
data <- get(load("E:/SHMC_1805/Lab data/Sunlab/bioinformatics/rawdata/scRNA_harmony_UMAP.Rdata"))
dataL <- list()
dataL <- SplitObject(data, split.by = "orig.ident") 
#选择其中一个样本开始分析，节约运存
for(i in 1:10){
pbmc <- dataL[[i]]
pbmc
## 寻找最优pK值
sweep.res.list <- paramSweep_v3(pbmc, PCs = 1:30, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)  
bcmvn <- find.pK(sweep.stats)
pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric()

## 排除不能检出的同源doublets，优化期望的doublets数量
DoubletRate = 0.078                     # 10000细胞对应的doublets rate是7.8%
homotypic.prop <- modelHomotypic(pbmc$seurat_clusters)   # 最好提供celltype
nExp_poi <- round(DoubletRate*ncol(pbmc)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## 使用确定好的参数鉴定doublets
pbmc <- doubletFinder_v3(pbmc, PCs = 1:30, pN = 0.25, pK = pK_bcmvn, 
                         nExp = nExp_poi.adj, reuse.pANN = F, sct = F)


## 结果展示，分类结果在pbmc@meta.data中,保存数据
pbmc$doubFind_res = pbmc@meta.data %>% select(contains('DF.classifications'))
pbmc$doubFind_score = pbmc@meta.data %>% select(contains('pANN'))
table(pbmc$DF.classifications_0.25_0.18_510)
scRNAlist <- list()
scRNAlist[[i]] <- pbmc@meta.data$doubFind_res %in% "singlet"}

#统计结果以及不同cluster的doubletscore
colnames(pbmc@meta.data)[ncol(pbmc@meta.data)]="DoubletFinder"
DF_df <- pbmc@meta.data[,c("seurat_clusters","DoubletFinder")]
write.table(DF_df, file = "DoubletFinder_result.txt", quote = FALSE, sep = '\t', row.names = FALSE, col.names = TRUE)
p1_doublet <- read.csv("DoubletFinder_result.txt",sep = '\t',stringsAsFactors = T)
p1_doublet$seurat_clusters <- as.character(p1_doublet$seurat_clusters)
p1_doublet$group <- "p1"
#系数分别为PN,PK,nEXP_Poi.Adj，展示umap图
DimPlot(pbmc,reduction = "umap",group.by ="DF.classifications_0.25_0.18_510")
#统计doubletscore
library(ggplot2)
ggplot(data=p1_doublet, aes(x=seurat_clusters,y=DoubletFinder,color=seurat_clusters))+
  geom_boxplot()+theme_classic()
ggsave('doubletscore_box_plot.png', width = 8.24, height = 2.84)