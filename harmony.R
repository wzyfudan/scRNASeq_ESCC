rm(list=ls())
setwd("E:/Lab data/Sunlab/bioinformatics/rawdata")
#if(!require(harmony))devtools::install_github("immunogenomics/harmony")
#devtools::install_github(repo = "mojaveazure/loomR")
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
BiocManager::install("harmony")
library(cowplot)
library(dplyr)
library(reshape2)
library(loomR)
library(patchwork)
library(Seurat)
library(tidyverse)
Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
package_version(harmony)

#整合数据利用list,paste命令构建路径变量dir(此处可跳过)
a = list.files("E:/SHMC_1805/Lab data/Sunlab/ESO scTCR-BCR/Rawdata/GEX/")                                                       
dir = paste("E:/SHMC_1805/Lab data/Sunlab/ESO scTCR-BCR/Rawdata/GEX/",a,sep="")
sample_name = c("P1", "P2","P3","P4","P5","P6","P7","P8","P9","P10")
#标准化路径处理，读取和过滤
scRNAlist <- list()
for(i in 1:length(dir)){
  counts <- Read10X(data.dir = dir[i])
  scRNAlist[[i]] <- CreateSeuratObject(counts, project=sample_name[i], min.cells=3, min.features = 200)
  scRNAlist[[i]][["percent.mt"]] <- PercentageFeatureSet(scRNAlist[[i]], pattern = "^MT-")
  scRNAlist[[i]] <- subset(scRNAlist[[i]], subset = nCount_RNA > 500 & nCount_RNA < 50000 & nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 15) 
}   
#多样本粗略merge
names(scRNAlist) <- sample_name

#此处可以导入经过doubletfinder处理过的数据
scRNAlist <- readRDS("scRNAlist_singletonly.rds")
scRNA <- merge(scRNAlist[[1]], y=c(scRNAlist[[2]], scRNAlist[[3]], 
                                   scRNAlist[[4]], scRNAlist[[5]], scRNAlist[[6]], scRNAlist[[7]], 
                                   scRNAlist[[8]], scRNAlist[[9]], scRNAlist[[10]]))
table(scRNA@meta.data$orig.ident)
#绘制质控图
plot.features = c("nFeature_RNA", "nCount_RNA", "percent.mt")
group = "orig.ident"
plots=list()
for(i in seq_along(plot.features)){
  plots[[i]] = VlnPlot(scRNA, group.by=group, pt.size = 0, features = plot.features[i]) + NoLegend()
}
plot = violin <- wrap_plots(plots = plots, nrow=2)
plot
dev.off()
#保存预处理数据
save(scRNA,file = 'scRNA.Rdata')

##==harmony整合多样本==##
rm(list=ls())
library(Rcpp)
library(harmony)
scRNA_harmony <- get(load("2023-1-6 T cell/CD8T_cell_harmony.Rdata"))
scRNA_harmony <- readRDS("B_cell.Rds")
scRNA_harmony <- readRDS("Monocytes_harmony.Rds")
scRNA_harmony <- readRDS("Macro.Rds")
scRNA_harmony <- readRDS("CD8T_cell.Rds")
scRNA_harmony <- scRNA_harmony[,scRNA_harmony@meta.data$seurat_clusters %in% c(0)]

##标准化+PCA降维
scRNA_harmony <- NormalizeData(scRNA_harmony) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose=FALSE)
scRNA_harmony <- scRNA_harmony %>% RunHarmony("orig.ident", plot_convergence = T)
#降维聚类，注意保持dims一致，调整res由大到小
scRNA_harmony <- RunUMAP(scRNA_harmony, reduction = "harmony", dims = 1:20)
scRNA_harmony <- FindNeighbors(scRNA_harmony, reduction = "harmony", dims = 1:20)
library(clustree) #不同 resolution 的树状图展示
{scRNA_harmony <- FindClusters(scRNA_harmony, resolution = 0.1)
scRNA_harmony <- FindClusters(scRNA_harmony, resolution = 0.2)
scRNA_harmony <- FindClusters(scRNA_harmony, resolution = 0.4)
scRNA_harmony <- FindClusters(scRNA_harmony, resolution = 0.6)
scRNA_harmony <- FindClusters(scRNA_harmony, resolution = 0.8)
scRNA_harmony <- FindClusters(scRNA_harmony, resolution = 1)
scRNA_harmony <- FindClusters(scRNA_harmony, resolution = 1.2)
clustree(scRNA_harmony)}
scRNA_harmony <- FindClusters(object = scRNA_harmony, verbose = T, resolution = 0.6)
save(scRNA_harmony, file = "scRNA_harmony_UMAP.Rdata")
saveRDS(scRNA_harmony, file = "CXCL13 Tex.Rds")
##作图
#group_by_cluster
plot1 = DimPlot(scRNA_harmony, reduction = "umap", label=T) 
#group_by_sample
plot2 = DimPlot(scRNA_harmony, reduction = "umap", group.by='group') 
plot1|plot2
library(ggplot2)
ggsave(filename = 'T cell_harmanoy.PDF', width = 10, height = 6)
gc()
#统计
Cluster_prop <- table(scRNA_harmony@meta.data[,c("seurat_clusters","group")]) %>% as.data.frame()
ggplot(Cluster_prop,aes(x=seurat_clusters,y=Freq,
                   fill=group))+  #根据不同Type填充颜色
  geom_bar(stat = 'identity')+
  theme_bw(base_size = 18)+ 
  theme(axis.text = element_text(colour = 'black'))

exhaust_markers <- list(c("HAVCR2", "PDCD1", "CTLA4","TIGIT","TOX2"))
cytotoxic_markers <- list(c("PRF1","NKG7","GNLY","GZMB","GZMA","TYROBP"))
eff_markers <- list(c("PRF1","IFNG","CCL4","GZMK"))
pbmc <- AddModuleScore(object = scRNA_harmony, features =exhaust_markers, name = "cell")
pbmc <- AddModuleScore(object = scRNA_harmony, features =c("GZMB"), name = "cell")
VlnPlot(pbmc,features = 'cell1',group.by = "seurat_clusters")+ylab("GZMK")+ggtitle("")


#统计
meta = scRNA_harmony@meta.data
NR = rbind(meta[which(meta$orig.ident == 'P1'),],meta[which(meta$orig.ident == 'P5'),],meta[which(meta$orig.ident == 'P7'),],meta[which(meta$orig.ident == 'P9'),],meta[which(meta$orig.ident == 'P10'),])
R = rbind(meta[which(meta$orig.ident == 'P2'),],meta[which(meta$orig.ident == 'P3'),],meta[which(meta$orig.ident == 'P4'),],meta[which(meta$orig.ident == 'P6'),],meta[which(meta$orig.ident == 'P8'),])
ids1 = table(NR[,c("seurat_clusters","orig.ident")]) %>% as.data.frame() 
ids2 = table(R[,c("seurat_clusters","orig.ident")]) %>% as.data.frame()
ids1$group = 'NR'
ids2$group = 'R'
df = rbind(ids1,ids2)
#去除值为0的行
df$Freq[df$Freq == 0]<-NA
df = na.omit(df)
# 使用频数计算比例
df$prop = 0
ESOtest1 <- readRDS("E:/Lab data/Sunlab/bioinformatics/rawdata/scRNA_harmony_UMAP.Rds")
meta = ESOtest1@meta.data
ids = table(meta$orig.ident) %>% as.data.frame()
for(i in 1:dim(ids)[1]){
  prop.temp = df[which(df$orig.ident %in% ids[i,1]),'Freq'] / ids[i,2]
  df[which(df$orig.ident %in% ids[i,1]),'prop'] = prop.temp
}
#df <- df %>%
#  group_by(orig.ident) %>%
#  mutate(total = sum(Freq),
#         prop = Freq / total)

#
library(rstatix)
df_pvals <- df %>%
  group_by(seurat_clusters) %>% #横坐标
  rstatix::wilcox_test(  #选择t_test或者wilcox或者kruskal_test()
    prop ~ group,  #纵坐标[重点更改]和分组
    p.adjust.method = "BH",  #还可以选择bonferroni
    ref.group = "NR") %>% 
  rstatix::add_significance("p") %>% #转化为
  rstatix::add_x_position(x = "seurat_clusters", dodge = 0.9)#[重点更改]确定p值的空间位置，否则会直接以分组显示或是p不显示
#write.csv(df_pvals,"df_pvals.csv")
#df_pvals <- read.csv("df_pvals.csv",row.names = 1)
#改动x,y的标签
#点图
library(ggprism)
library(ggpubr)

#重命名seurat_clusters
#df$seurat_clusters = factor(df$seurat_clusters, labels = c("CXCL13 Tex", "IL7R Tm", "TNFSF9 Tac",  "GZMK Teff", "CD74 Teff",  
#                                                           "ITGAE TRM", "FCGR3A Tc", "NK/NKT", "CX3CR1 Teff","ISG15 Tex")) 

#Boxplot,P值可以选择"p.signif"或者p.adj用*表示
#df_pvals <- subset(df_pvals, p < 0.06)

library(ggplot2)
ggplot(df, aes(x=seurat_clusters,y=prop))+geom_boxplot(aes(fill = group))+ 
  theme_bw()+ stat_pvalue_manual(df_pvals,label = "p",tip.length = 0.01,y.position = 0.06)
ggsave("CD4 T_wilcox.PDF",width = 10,height = 6)
