
rm(list=ls())
setwd("E:/Lab data/Sunlab/bioinformatics/rawdata")
#if (!require("BiocManager", quietly = TRUE))
  #install.packages("BiocManager")
#BiocManager::install("harmony")
library(dplyr)
library(Seurat)
library(patchwork)
library(tidyverse)
library(ggplot2)
#整合数据利用list,paste命令构建路径变量dir
a = list.files("E:/Lab data/Sunlab/ESO scTCR-BCR/Rawdata/GEX/")                                                       
dir = paste("E:/Lab data/Sunlab/ESO scTCR-BCR/Rawdata/GEX/",a,sep="")
sample_name = c("P1", "P2","P3","P4","P5","P6","P7","P8","P9","P10")
#标准化路径处理，读取和过滤
scRNAlist <- list()
for(i in 1:length(dir)){
  counts <- Read10X(data.dir = dir[i])
  scRNAlist[[i]] <- CreateSeuratObject(counts, project=sample_name[i], min.cells=3, min.features = 200)
  scRNAlist[[i]][["percent.mt"]] <- PercentageFeatureSet(scRNAlist[[i]], pattern = "^MT-")}
scRNA <- merge(scRNAlist[[1]], scRNAlist[2:length(scRNAlist)]) #merge自动添加barcode的样本下标，区分批次
plots1 <- FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "percent.mt")
plots2 <- FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plots1+plots2
table(scRNA$orig.ident)
for(i in 1:length(dir)){
  scRNAlist[[i]] <- subset(scRNAlist[[i]], subset = nCount_RNA > 500 & nCount_RNA < 50000 & nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 15) 
}
#查看细胞数
names(scRNAlist) <- sample_name
scRNA <- merge(scRNAlist[[1]], scRNAlist[2:length(scRNAlist)])
table(scRNA$orig.ident)
saveRDS(scRNAlist,file = 'scRNA_orig.Rds')
#绘制质控图
plot.features = c("nFeature_RNA", "nCount_RNA", "percent.mt")
group = "orig.ident"
plots=list()
for(i in seq_along(plot.features)){
  plots[[i]] = VlnPlot(scRNA, group.by=group, pt.size = 0, features = plot.features[i]) + NoLegend()
}
plot = violin <- wrap_plots(plots = plots, nrow=2)
plot
