#https://www.jianshu.com/p/75c6d55a63bb?from=timeline
#https://blog.csdn.net/qq_45478665/article/details/119891954
#http://events.jianshu.io/p/159df5c4ff21

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


#integrte整合分析计算量较大:http://events.jianshu.io/p/159df5c4ff21
rm(list=ls())
scRNAlist <- readRDS("E:/Lab data/Sunlab/bioinformatics/rawdata/scRNA_orig.Rds")
#删除整理内存，并升级运存，设置future并行http://www.360doc.com/content/21/0715/15/76149697_986667695.shtml
library(future)
options(future.globals.maxSize = 50 * 1024^3)
#SCTransform标准化
scRNAlist <- parallel::mclapply(scRNAlist, FUN=function(x) SCTransform(x), mc.cores=1)
#寻找前3000的可变基因
scRNA.features <- SelectIntegrationFeatures(object.list = scRNAlist, nfeatures = 3000)
scRNAlist <- PrepSCTIntegration(object.list = scRNAlist, anchor.features = scRNA.features, 
                                verbose = FALSE)
#寻找anchors，计算量巨大
scRNA.anchors <- FindIntegrationAnchors(object.list = scRNAlist, normalization.method = "SCT", 
                                          anchor.features = scRNA.features,dims = 1:30)
scRNA.integrated <- IntegrateData(anchorset = scRNA.anchors, normalization.method = "SCT",dims = 1:30)
# 切换至IntegrateData整合分析
DefaultAssay(ESOtest.integrated) <- "integrated"
# 运行标准分析
scRNA.integrated <- ScaleData(scRNA.integrated, verbose = FALSE)
scRNA.integrated <- RunPCA(scRNA.integrated, npcs = 50, verbose = FALSE)
#elbow确定维度，选择10个主成分代表数据集
ElbowPlot(scRNA, ndims=50)
pc.num=1:10
scRNA <- RunUMAP(scRNA, dims=pc.num)
PLot1 <- DimPlot(scRNA.integrated, reduction = "tsne", group.by = "orig.ident")
PLot2 <- DimPlot(scRNA, group.by = "orig.ident", split.by = "orig.ident", ncol =4)
dev.off()
#保存聚类结果
save(scRNA,file="scRNA_UMAP.rda")



#单样本步骤详解
# 使用Read10X函数读取矩阵数据，得到一个稀疏矩阵
counts <- Read10X(data.dir = dir)
ESOtest = CreateSeuratObject(counts, project = "ESOtest", min.cells = 3, min.features = 200)
#查看每个样本的细胞数
table(ESOtest@meta.data$orig.ident) 
#使用从MT-作为线粒体基因集
ESOtest[["percent.mt"]] <- PercentageFeatureSet(ESOtest, pattern = "^MT-")
#qc指标可视化
VlnPlot(ESOtest, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#过滤具有超过2500或少于200个基因计数的细胞，同时过滤掉线粒体比例超过5%的细胞
ESOtest <- subset(ESOtest, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
#通过总表达值对每个细胞的基因表达值归一化，并将其乘以缩放因子
ESOtest <- NormalizeData(ESOtest, normalization.method = "LogNormalize", scale.factor = 10000)
#ESOtest[["RNA"]]@data[c("CD3D", "TCL1A", "MS4A1"), 1:10]
#高变异基因选择,计算数据集中表现出高细胞间差异的基因子集,这里选择2000个
ESOtest<- FindVariableFeatures(ESOtest, selection.method = "vst", nfeatures = 2000)
# 前10的高变异基因
top10 <- head(VariableFeatures(ESOtest), 10)
# 将前10的高变异基因在图中标注出来
plot1 <- VariableFeaturePlot(ESOtest)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2
#数据缩放
all.genes <- rownames(ESOtest)
ESOtest <- ScaleData(ESOtest, features = all.genes)
#PCA降维分析
ESOtest <- RunPCA(ESOtest, features = VariableFeatures(object = ESOtest))
#选择前5个维度进行查看
print(ESOtest[["pca"]], dims = 1:5, nfeatures = 5)
#作PCA图
DimPlot(ESOtest, reduction = "pca")
DimHeatmap(ESOtest, dims = 1:5, cells = 500, balanced = TRUE)
#维度确认，选择10个主成分代表数据集，这里选择10个后续作图pc.num都是1：10
#将数据的1%打乱重新运行pca，构建特征得分的“零分布”，这个过程重复100次
ESOtest <- JackStraw(ESOtest, num.replicate = 100)
ESOtest <- ScoreJackStraw(ESOtest, dims = 1:20)
plot1<-JackStrawPlot(ESOtest, dims = 1:15)
plot2<-ElbowPlot(ESOtest)
plot1
#以上为质控部分
#1.QC阶段（关键步骤），nUMI（nCount_RNA),nGene(nFeature_RNA),percent.ribo线粒体比例的筛选
#2.多样本merge和integrate阶段，影响整合强度等参数的设置以及整合方法的选择（seurat自带的'FindIntegrationAnchors还是Harmony等方法）
#3.ScaleData参数设置，输入的基因''features''(可变基因或全部基因），需要回归的因素''vars.to.regress''（UMI和percent.ribo）都会影响后期的UMAP结果
#4.PCA数量的选择，选择合适的npca数量（elbow图的拐角处对应的参数，尝试多个npca作UMAP图找到最佳分群的npca）作为FindNeighbors的输入
#5.FindClusters阶段，resolution参数（分辨率）的选择决定分群的数量，可以通过''clustree''包查看不同res的分群情况
#查看之前输出数据结果
table(ESOtest@meta.data$orig.ident)

#细胞聚类
#计算最邻近距离
ESOtest <- FindNeighbors(ESOtest, dims = 1:10)
#聚类，包含设置下游聚类的“间隔尺度”的分辨率参数resolution ，增加值会产生更多聚类Cluster。
ESOtest <- FindClusters(ESOtest, resolution = 0.5)
# 查看前5个细胞的聚类id
head(Idents(ESOtest), 5)
table(ESOtest@meta.data$seurat_clusters)
metadata <- ESOtest@meta.data
cell_cluster <- data.frame(cell_ID=rownames(metadata), cluster_ID=metadata$seurat_clusters)
#非线性降维
#tSNE
ESOtest = RunTSNE(ESOtest, dims = 1:10)
embed_tsne <- Embeddings(ESOtest, 'tsne')   #提取tsne图坐标
#group_by_cluster
plot1 = DimPlot(ESOtest, reduction = "tsne", label=T) 
#group_by_sample
plot2 = DimPlot(ESOtest, reduction = "tsne", group.by='orig.ident') 
#combinate
plotc <- plot1+plot2
plotc
#UMAP
ESOtest <- RunUMAP(ESOtest, dims = 1:10)
embed_umap <- Embeddings(ESOtest, 'umap')   #提取umap图坐标
#group_by_cluster
plot3 = DimPlot(ESOtest, reduction = "umap", label=T) 
#group_by_sample
plot4 = DimPlot(ESOtest, reduction = "umap", group.by='orig.ident')
#combinate
plot3
