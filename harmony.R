
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
