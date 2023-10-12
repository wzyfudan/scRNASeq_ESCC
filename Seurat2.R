rm(list=ls())
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
setwd("E:/Lab data/Sunlab/bioinformatics/rawdata")
ESOtest <- readRDS("E:/Lab data/Sunlab/bioinformatics/rawdata/B_cell_harmony.Rds")
ESOtest <- readRDS("E:/Lab data/Sunlab/bioinformatics/rawdata/scRNA_harmony_UMAP.Rds")
ESOtest <- readRDS("E:/Lab data/Sunlab/bioinformatics/rawdata/CD4T_cell.Rds")
ESOtest <- readRDS("E:/Lab data/Sunlab/bioinformatics/rawdata/CXCL13 Tex.Rds")
table(ESOtest@meta.data$seurat_clusters)
#寻找cluster之间差异基因
#找到每个cluster相比于其余cluster的markgene，自定义logfc；test.use是差异marker的检验方式，默认为"wilcox",可选"MAST"，"DESeq2"等
ESOtest.markers <- FindAllMarkers(ESOtest, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
#table使用交叉分类的因素建立一个列联表,计数在每个组合的因素水平。
table(ESOtest.markers$cluster)
saveRDS(ESOtest.markers,file = 'ESO_markers.Rds')
ESOtest.markers <- readRDS("E:/Lab data/Sunlab/bioinformatics/rawdata/2022-12-2 doublet+Findmarkers/ESO_markers.Rds")
ESOtest.markers <- readRDS("E:/Lab data/Sunlab/bioinformatics/rawdata/2023-1-6 T cell/T_markers.Rds")
#绘制每个cluster按logfc排名前10的markgene的热图
top10 <- ESOtest.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.csv(top10,"T_top10.csv")
write.csv(ESOtest.markers,"ALL_markers.csv")
#data给scale.data赋值
DoHeatmap(object = ESOtest, features = top10$gene, slot = "scale.data",label = F, group.bar = T, group.bar.height = 0.03)
ggsave('Doheatmap_markers.PDF', width = 12, height = 10)
dev.off()
#查看特定基因在聚类图中的表达量featureplot分布情况,TBSAB1-Mast,GNLY-NK/NKT,LYZ-Myeloid，IL3RA-pDC
FeaturePlot(ESOtest, features = c("CD8A","GZMK","GZMB","PDCD1","CTLA4","CXCL13","ITGAE","FCGR3A","CD27"))
FeaturePlot(ESOtest, features = c("CD4", "IL7R","FOXP3","CXCL13","IFNG","HLA-DRA","PDCD1","CTLA4","TNFRSF9"))
FeaturePlot(ESOtest, features = c("CD79A","MS4A1","JCHAIN","MZB1","IGHD","CD27","BCL6","TCL1A","FCRL5"))
FeaturePlot(ESOtest, features = c("CD68","CD163","CD14","FCGR1A","S100A8","CSF3R"))
ggsave('Featureplot_CD4.PNG', width = 12, height = 9)
dev.off()

#细胞手动注释
new.cluster.ids <- c("CD8 T cell", "B cell", "CD4 T cell",  "Myeloid cell", "CD8+NKT-like cell",  
                     "Plasma B cell", "Cycling T cell", "Mast cell", "Native B","Fibroblast","Epithelial cell","Fibroblast")
new.cluster.ids <- c("CXCL13 Tex", "IL7R Tm", "TNFSF9 Tcm",  "GZMK Tem", "CD74 Tem",  
                     "ITGAE Trm", "FCGR3A NK-T", "TRDV1 gdT", "CX3CR1 Temra","ISG15 Tex")
new.cluster.ids <- c("IL7R Tm", "LAIR2 Treg", "IFNG Th1/Tfh",  "RTKN2 Treg", "CXCL13 Th1/Tfh", "KLF2 Tm",  
                     "ISG15 Treg", "GZMA Tem", "TNFRSF9 Treg", "TNFRSF4 Treg")
new.cluster.ids <- c("TNFRSF13B Bmem","FCER2 Bn","IGHG PC","B/T mix","FCRL4 Bmem","NR4A1 Bmem","IGHG PC","TCL1A BGC","STMN1 BGC","IGHG PC","IGHG PC","IGHA PC")
new.cluster.ids <- c("S100A8 Neu","CD8 T","TPSAB1 Mast","CD14 Mono","PLPP3 Neu","C1QC Macro","CD1c cDC2","CCL3L1 Neu","ISG15 Neu","LILRA4 pDC","LAMP3 DC3","CD79A B","CLEC9A cDC1")
#new.cluster.ids <- c("CXCR4 Teff","TNFAIP3 T","ZNF683 Trm","XCL1 T","NKG7 Tc","CXCL13 T","TIGIT Tex","GZMK Tem")
names(new.cluster.ids) <- levels(ESOtest)
ESOtest <- RenameIdents(ESOtest, new.cluster.ids)
ESOtest@meta.data$new.cluster.ids <- Idents(ESOtest)#并将信息添加至metadata
{ESOtest$group <- "R" 
ESOtest@meta.data$group[ESOtest@meta.data$orig.ident %in% c("P1","P5","P7","P9","P10")] = "NR"}
{ESOtest$TRG <- "TRG I" 
ESOtest@meta.data$TRG[ESOtest@meta.data$orig.ident %in% c("P2","P5","P6")] = "TRG II"
ESOtest@meta.data$TRG[ESOtest@meta.data$orig.ident %in% c("P1","P7","P9","P10")] = "TRG III"}
#作图并显示各组分B细胞比例***
df = data.frame(clu=names(table(ESOtest$seurat_clusters)),
                per=sprintf("%1.2f%%", 100*table(ESOtest$seurat_clusters)/length(ESOtest$seurat_clusters)))
ESOtest$per = df[match(ESOtest$seurat_clusters,df$clu),2]
ESOtest$new = paste0(ESOtest$seurat_clusters,":",ESOtest$new.cluster.ids,"(",ESOtest$per,")")
DimPlot(ESOtest,reduction = "umap",group.by = "new.cluster.ids",label = T,repel = T)+guides(color = FALSE) +ggtitle("Myeloid cell")
DimPlot(ESOtest,reduction = "umap",group.by = "TRG",label = F)+
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank(),
        panel.border = element_blank())
ggsave('Mye_annotation.PDF', width = 9, height = 6)
ggsave('Mono.Tiff', width = 10, height = 6)
saveRDS(ESOtest,file = "B_cell_harmony.Rds")

#完成按照细胞类别的表达热图heatmap
ESOtest.markers <- readRDS("2022-12-2 doublet+Findmarkers/ESO_markers.Rds")
top10 <- ESOtest.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(object = ESOtest, features = top10$gene, label = F)
ggsave('heatmap.PDF', width = 10, height = 10)
#高级作图，火山图
install.packages('devtools')
devtools::install_github("jokergoo/ComplexHeatmap")
devtools::install_github("sajuukLyu/ggunchull", type = "source")#可以有效的圈起来需要的群
devtools::install_github("junjunlab/jjAnno") #轻松给堆积条形图添加注释
devtools::install_github('junjunlab/scRNAtoolVis')
library(scRNAtoolVis)
#ESOtest.markers$cluster <- factor(ESOtest.markers$cluster, levels = c("0", "1","2","3","4","5","6","7","8","9","10","11"), labels = c("CD8 T cell", "B cell", "CD4 T cell",  "Myeloid cell", "CD8+NKT-like cell",  
                                                                                                            #  "Plasma B cell", "Cycling T cell", "Mast cell", "Native B","Fibroblast","Epithelial cell","Fibroblast"))
jjVolcano(diffData = ESOtest.markers,
          log2FC.cutoff = 0.5, 
          size  = 3.5, #设置点的大小
          #col.type = "adjustP",
          fontface = 'italic', #设置字体形式
          #aesCol = c('purple','orange'), #设置点的颜色
          tile.col = corrplot::COL2('RdBu', 15)[2:13], #设置cluster的颜色,数量设置对应
          #col.type = "adjustP", #设置矫正方式
          topGeneN = 3, #设置展示topN的基因
          legend.position = c(0.8,1),
)
ggsave('jjVolcano.PDF', width = 14, height = 6)
dev.off()


#挑选感兴趣的子集(2，8条）——https://cloud.tencent.com/developer/article/2036759
sce <- get(load("E:/Lab data/Sunlab/bioinformatics/rawdata/scRNA_harmony_UMAP.Rdata"))
table(sce@meta.data$seurat_clusters)
#提取子集, 需要知道逻辑值，坐标或名字
# seurat 取聚类后的0,2,4,6分群,无需重新聚类//subset需要重新聚类
T_sce <- sce[,sce@meta.data$seurat_clusters %in% c(0,2,4,6)]
B_sce <- sce[,sce@meta.data$seurat_clusters %in% c(1,5,8)]
Mono_sce <- sce[,sce@meta.data$seurat_clusters %in% c(3,7)]
CD8T_sce <- sce[,sce@meta.data$seurat_clusters %in% c(0,2,6,7)]
CD4T_sce <- sce[,sce@meta.data$seurat_clusters %in% c(1,3,4)]
Treg <- ESOtest[,ESOtest@meta.data$seurat_clusters %in% c(1,3,6,8,12)]
sce <- ESOtest[,ESOtest@meta.data$seurat_clusters %in% c(0,2,3,4,5,6,7,8,9,10,12)]
sce <- ESOtest[,ESOtest@meta.data$seurat_clusters %in% c(3,5)]
saveRDS(sce,"Macro.Rds")

DimPlot(sce, reduction = "umap", label = T,pt.size = 0.5)
#保存亚集
saveRDS(sce,file = "Mono_cell.Rds")
save(Treg,file = "Treg_cell.Rdata")
#或者使用名称提取如 cd4_sce <- sce[, Idents(sce) %in% c("Naive CD4 T" ,  "Memory CD4 T" )]
#取子集后就是新的项目，需要重新标准化normalization到FindVaribaleFeatures(harmony)


#差异分析，在同一cluster不同分组之间找差异基因
#注释分组信息
seurat@meta.data$new.cluster.ids <- Idents(seurat)
seurat$group <- "R" 
seurat@meta.data$group[seurat@meta.data$orig.ident %in% c("P1","P5","P7","P9","P10")] = "NR"
#差异分析,选定某个cluster
i=c("7")
subobj <- subset(ESOtest@meta.data,seurat_clusters %in% i)
scRNAsub <- subset(ESOtest, cells=row.names(subobj))
DefaultAssay(scRNAsub) <- 'RNA'
dge.celltype <- FindMarkers(scRNAsub, ident.1 = 'R',ident.2 = 'NR', 
                            group.by = 'group',logfc.threshold = 0.05,only.pos = F) #1是case,2是control,结果pct1和pct2分别表示两组表达细胞的占比,avg_logFC两组间平均logFC默认取log2
dge.celltype <- FindMarkers(ESOtest, ident.1 = 'R',ident.2 = 'NR', 
                            group.by = 'group',logfc.threshold = 0.05,only.pos = F) #总分析
write.csv(dge.celltype,"dge.celltype.csv")
#AverageExpression验证log2FC
av <-AverageExpression(scRNAsub, assays="RNA", group.by="group",verbose=F) 
av$RNA[c("TSC22D3","DUSP1","GNLY"),]

#哪些marker基因在不同的样本中“守恒”，相当于差异分析的反面
#BiocManager::install('multtest')
#install.packages('metap')
#dge.celltype <- FindConservedMarkers(ESOtest, ident.1 = 1,grouping.var = "group",logfc.threshold = 0.25,min.diff.pct = 0.25,min.pct = 0.25,only.pos = T)
#火山图可视化
#BiocManager::install('EnhancedVolcano')
library("EnhancedVolcano")
library("ggrepel")
dge.celltype <- read.csv("pathway/NR_R DEGs logFC1&pvalue0.05.csv",row.names = 1)
EnhancedVolcano (dge.celltype,
                 lab=rownames(dge.celltype),
                 x="log2FoldChange", #红色表示同时满足P值和表达量两个筛选条件的基因
                 y="pvalue", #Benjamini-Hochberg校正方法
                 pCutoff=0.05,
                 FCcutoff=3,
                 pointSize=2.0,
                 xlab = bquote("Responsive" ~ "--" ~ Log[2] ~ "fold change" ~ "--" ~"Non-Responsive"),     # 此行为默认x轴名
                 ylab = bquote(~-Log[10] ~ "p-value"), # 修改y轴名为斜体
                 title="R VS NR",
                 subtitle=bquote(italic("DEGs in bulk-seq")),
                 titleLabSize=15,
                 ylim = c(1,6),
                 selectLab = c("S100A7A","S100A7","CXCL9"))+ theme(axis.text.x = element_text(color="black", size=9),
                                         axis.text.y = element_text(color="black", size=9),
                                         plot.title = element_text(hjust = 0)) 
ggsave('deg_S100A7_Volcano.PDF', width = 8, height = 6)


##患者入组group设置绘图
df <- data.frame(group = c("P1", "P5", "P7", "P9", "P10","P2","P3","P4","P6","P8"), value = c(90, 40, 80, 65, 60,10,0,0,15,0))
my_colors <- c(rep("blue", 5), rep("red", 5))
df$group <- factor(df$group, levels = c("P1", "P5", "P7", "P9", "P10","P2","P3","P4","P6","P8"))
ggplot(df, aes(x = group, y = value, fill = my_colors)) +
  geom_bar(stat = "identity", width = 0.5) +
  xlab("non-respond" ~ "←" ~ "Patients" ~ "→" ~"respond") +
  ylab("p Tumor Residual Ratio") +geom_hline(yintercept = c(1,50), linetype = "dashed")+
  geom_text(x=10,y=-0.5,label = "TRG I")+geom_text(x=10,y=25,label = "TRG II")+geom_text(x=10,y=75,label = "TRG III")+
  theme(legend.position = "none")
ggsave("Patient group.pdf",width = 8,height = 6)
##细胞比例按患者设置绘图
df <- table(ESOtest@meta.data[,c("new.cluster.ids","orig.ident")]) %>% as.data.frame()
library(dplyr)
df <- df %>%
  group_by(orig.ident) %>%
  mutate(Prop = Freq / sum(Freq))
p <- ggplot(df, aes(x = orig.ident, y = Prop, fill = new.cluster.ids)) +
  geom_bar(stat = "identity",  position = "stack") #堆叠区别dodge
p
ggsave("Patient cluster prop_stack.PDF",width = 10,height = 6)
