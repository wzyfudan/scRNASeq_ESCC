#https://www.jianshu.com/p/2805740c1624
#https://ncborcherding.github.io/vignettes/vignette.html#1_Introduction
#https://blog.csdn.net/weixin_42382703/article/details/112511753
rm(list=ls())
setwd("E:/Lab data/Sunlab/bioinformatics/rawdata")
# Most up-to-date version
#devtools::install_github("ncborcherding/scRepertoire@dev")
#BiocManager::install("scRepertoire",force = TRUE)
#packageVersion("scRepertoire")
# 加载所需的R包
#BiocManager::install("GenomeInfoDbData")
library("GenomeInfoDbData")
library(tidyverse)
suppressMessages(library("scRepertoire"))
packageVersion("scRepertoire") #‘1.10.1’
suppressMessages(library("Seurat"))
library(ggplot2)
library(ggpubr)

# 读取不同的样本contig文件
S1 <- read.csv("E:/Lab data/Sunlab/data/ESO scTCR-BCR/Rawdata/TCR/TCR/P1_filtered_contig_annotations.csv")
S2 <- read.csv("E:/Lab data/Sunlab/data/ESO scTCR-BCR/Rawdata/TCR/TCR/P2_filtered_contig_annotations.csv")
S3 <- read.csv("E:/Lab data/Sunlab/data/ESO scTCR-BCR/Rawdata/TCR/TCR/P3_filtered_contig_annotations.csv")
S4 <- read.csv("E:/Lab data/Sunlab/data/ESO scTCR-BCR/Rawdata/TCR/TCR/P4_filtered_contig_annotations.csv")
S5 <- read.csv("E:/Lab data/Sunlab/data/ESO scTCR-BCR/Rawdata/TCR/TCR/P5_filtered_contig_annotations.csv")
S6 <- read.csv("E:/Lab data/Sunlab/data/ESO scTCR-BCR/Rawdata/TCR/TCR/P6_filtered_contig_annotations.csv")
S7 <- read.csv("E:/Lab data/Sunlab/data/ESO scTCR-BCR/Rawdata/TCR/TCR/P7_filtered_contig_annotations.csv")
S8 <- read.csv("E:/Lab data/Sunlab/data/ESO scTCR-BCR/Rawdata/TCR/TCR/P8_filtered_contig_annotations.csv")
S9 <- read.csv("E:/Lab data/Sunlab/data/ESO scTCR-BCR/Rawdata/TCR/TCR/P9_filtered_contig_annotations.csv")
S10 <- read.csv("E:/Lab data/Sunlab/data/ESO scTCR-BCR/Rawdata/TCR/TCR/P10_filtered_contig_annotations.csv")
# 构建重叠群列表
contig_list <- list(S1,S2,S3,S4,S5,S6,S7,S8,S9,S10)
# 查看数据集中contig的注释信息
head(contig_list[[1]])
#合并VDJ重叠群，构建TCR克隆型,将细胞的barcode与clonotype直接对应起来,同时增加ID以防重复
combined <- combineTCR(contig_list, 
                       samples = c("P1","P2","P3","P4","P5","P6","P7","P8","P9","P10"), 
                       ID = c("NR","R","R","R","NR","R","NR","R","NR","NR"), 
                       removeNA = TRUE,
                       removeMulti = TRUE) #去除TCRA和TCRB多链或缺失者
#生成的结果有：核苷酸序列 (CTnt)、氨基酸序列 (CTaa)、VDJC 基因序列 (CTgene) 
#去除barcode前缀注释信息
for (i in 1:10) {
  combined[[i]] <- stripBarcode(combined[[i]], column = 1, connector = "_", num_connects = 3)
  }
#统一与seurat转录组一致的barcode，方便后续比较
df <- data.frame()
for (i in 1:10) {
  combined[[i]]$barcode <- gsub('-.*',paste0("-1_", i),combined[[i]]$barcode)
  df <- rbind(combined[[i]],df)
}
write.csv(df,"combinedTCR.csv")
df <- read.csv("E:/Lab data/Sunlab/bioinformatics/rawdata/T cell/2023-1-17 TCR_scRepertoire/combinedTCR.csv")
combined <- split(df,df$sample)
#去重
#df2 <- df[!duplicated(df[,1]),]
#combined <- split(df2,df2$sample)



#后续分析调用cloneCall参数：
#"gene"：使用包含TCR/Ig的VDJC基因
#"nt"：使用CDR3区域的核苷酸序列; "aa"：使用CDR3区域的氨基酸序列
#"gene+nt"：使用包含TCR/Ig的VDJC基因+CDR3区域的核苷酸序列
#PS:gene方法将是最敏感的，而使用nt或aa则是适度敏感的，并且对克隆型最具特异性的则是gene+nt。
#可视化高表达的clonetype的基因组成
vizGenes(combined[c(1,5,7,9,10)], gene = "V", 
         chain = "TRB", 
         plot = "bar", 
         order = "variance", 
         scale = TRUE)
vizGenes(combined[c(2,3,4,6,8)], gene = "V", 
         chain = "TRB", 
         plot = "bar", 
         order = "variance", 
         scale = TRUE)
vizGenes(combined[c(1,5,7,9,10)], 
         gene = "V", 
         chain = "TRB", 
         y.axis = "J", 
         plot = "heatmap", 
         scale = TRUE, 
         order = "gene")
##根据患者信息分析
#可视化,克隆size
quantContig(combined, cloneCall="gene+nt", scale = TRUE)
quantContig(combined, cloneCall="gene", group = "ID", scale = TRUE)
#设置输出clonaltype size表格结果
quantContig_output <- quantContig(combined, cloneCall="gene+nt", 
                                  scale = TRUE, exportTable = TRUE)
quantContig_output
#Clonal Proportion克隆型比例，1:10表示每个样本中的前10个克隆型。
clonalHomeostasis(combined, cloneCall = "aa")
clonalProportion(combined, cloneCall = "gene") 
#Overlap Analysis重叠分析
clonalOverlap(combined, cloneCall = "gene+nt", method = "morisita")
clonalOverlap(combined, cloneCall="aa", method="overlap")
#比较两者
compareClonotypes(combined, numbers = 10, 
                  samples = c("P8","P10"), 
                  cloneCall="aa", graph = "alluvial")
#Diversity Analysis多样性分析
clonalDiversity(combined, cloneCall = "gene", group = "samples", 
                n.boots = 100)
dev.off()

#---------------------------------------------------------------------
#与seurat联合分析
# 加载scRNA-seq数据
seurat <- readRDS("E:/Lab data/Sunlab/bioinformatics/rawdata/T_cell_harmony.Rds")
seurat <- readRDS("E:/Lab data/Sunlab/bioinformatics/rawdata/CD4T_cell.Rds")
seurat <- readRDS("E:/Lab data/Sunlab/bioinformatics/rawdata/CXCL13 Tex.Rds")
# 查看细胞降维聚类信息
DimPlot(seurat, group.by = "seurat_clusters",label = T) + NoLegend()
#combine整合,需要barcode去重
seurat <- combineExpression(combined, seurat,
                            cloneCall="gene", 
                            group.by = "none", 
                            proportion = FALSE, 
                            cloneTypes=c(Single=1, Small=5, Medium=20, Large=100, Hyperexpanded=1000))
#seurat1 <- seurat[,seurat$group %in% "NR"]
#seurat2 <- seurat[,seurat$group %in% "R"]
#统计NA即匹配不成功的数量
sum(is.na(seurat@meta.data$CTaa))

##统计分析
#使用内置的STARTrac——Zhang.2018
StartracDiversity(seurat, type = "seurat_clusters", sample = "orig.ident", by = "overall")
ggsave('CD4_STARTRAC.Tiff', width = 12, height = 9)

#计算repertoire相关参数，从seurat获得合并信息
combined2 <- expression2List(seurat, split.by = "new.cluster.ids") #cluster之间
combined2 <- expression2List(seurat, split.by = "seurat_clusters")
combined2 <- expression2List(seurat, split.by = "orig.ident")
combined2 <- expression2List(seurat, split.by = "group") #组与组之间
combined2 <- expression2List(seurat, split.by = "TRG")
#直方图展示各cluster中的扩增克隆型比例clonal size
occupiedscRepertoire(seurat, x.axis = "new.cluster.ids")
occupiedscRepertoire(seurat, x.axis = "orig.ident")+theme(axis.text.x = element_text(face = "bold", size = 12, hjust = 1),axis.text.y = element_text(face = "bold", size = 10, hjust = 1))
occupiedscRepertoire(seurat, x.axis = "group")
table <- occupiedscRepertoire(seurat, x.axis = "orig.ident",exportTable = TRUE) #叠加图
write.csv(table,"CD4_clonetype.csv")
ggsave('CD8_clonal patient.PDF', width = 10, height = 8)
#丰度Abundance
clonalHomeostasis(combined2, cloneCall = "nt") #组成比例
#Unique clonaltype 含义是Nr of clonotype/ Cell Number 
quantContig(combined2, cloneCall="gene+nt", scale = TRUE) 
sapply(unique(seurat@meta.data$orig.ident), function(x) length(unique(seurat@meta.data$CTaa[seurat@meta.data$orig.ident == x])))
#clonality & diversity(bootstrap估计,group.by是图例点，x.axis是互相比较的组)
clonalDiversity(combined2, cloneCall = "nt",group.by ="orig.ident",x.axis = "new.cluster.ids")
clonalDiversity(combined2, cloneCall = "nt",group.by ="orig.ident",x.axis = "TRG")+theme(axis.text.x = element_text(face = "bold", size = 12, angle = 45, hjust = 1))+
  stat_compare_means(comparisons = list(c("TRG I","TRG III"),c("TRG II","TRG III"))) + ylab("Clonal Diversity")+theme(axis.title.y = element_text(size = 12),axis.text.y = element_text(size = 11, hjust = 1))
ggsave('CD4_clonal diversity.PDF', width = 10, height = 8)

diversity_output <- clonalDiversity(combined2, cloneCall="nt",exportTable = TRUE,group.by ="new.cluster.ids",x.axis = "orig.ident")# exportTable =T输出表格计算P值
diversity_output <- clonalDiversity(combined2, cloneCall="nt",exportTable = TRUE,group.by ="orig.ident",x.axis = "group") #计算diversity
write.csv(diversity_output,"CD4_cluster TCR.csv")
#clonal sharing
clonalOverlap(combined2, cloneCall="aa", method="morisita")+
  theme(axis.text.x = element_text(face = "bold", size = 12, angle = 45, hjust = 1))+
  theme(axis.text.y = element_text(face = "bold",size = 12, hjust = 1))
ggsave('CD4T_clonal sharing.PDF', width = 7, height = 6)

##UMAP可视化-------------------------------------------------------------
#查看分布
slot(seurat, "meta.data")$cloneType <- factor(slot(seurat, "meta.data")$cloneType, 
                                              levels = c("Hyperexpanded (100 < X <= 1000)", 
                                                         "Large (20 < X <= 100)", 
                                                         "Medium (5 < X <= 20)", 
                                                         "Small (1 < X <= 5)", 
                                                         "Single (0 < X <= 1)", NA))
colorblind_vector <- colorRampPalette(c("#FF4B20", "#FFB433", 
                                        "#C6FDEC", "#7AC5FF", "#0348A6"))
DimPlot(seurat, group.by = "cloneType") +
  scale_color_manual(values = colorblind_vector(5), na.value="grey")

#高亮特征分布群，可选择 "CTaa" or "CTnt"，根据TCR CDR3序列
#{df_spe <- read.csv("mcpas_CD8.csv") #VDJbd,mcpas
#TRAB <- df_spe[,c("cdr3_aa1","cdr3_aa2")]}
{df_spe <- read.csv("2023-2-27 MANA score/GLIPH2_result.csv")
TRAB <- df_spe[,c("TcRa","TcRb")]}
TRAB <- TRAB[!apply(is.na(TRAB) | TRAB == "", 1, all),] #去除空白行
TRAB <- unite(TRAB, "CTaa",TcRa,TcRb, sep = "_", remove = FALSE) #合并TRA_TRB or cdr3_aa1,aa2
TRAB <- TRAB[1:1000,] %>% as.data.frame() #只选择前1000行
seurat <- highlightClonotypes(seurat, cloneCall= "aa", sequence = TRAB$CTaa) #报错重启
DimPlot(seurat, group.by = "highlight") + NoLegend() 
{CD8_spe <- table(seurat@meta.data[,c("new.cluster.ids","highlight")]) %>% as.data.frame()
aggregate(CD8_spe$Freq, by=list(type=CD8_spe$new.cluster.ids),sum)} #求特异性T在cluster分布数
# 按一定freq扩增频率，展示克隆分布等高线图
clonalOverlay(seurat, reduction = "umap", 
              freq.cutpoint = 30, bins = 10) + guides(color = "none")
#用箭头表示network interaction of clonotypes shared between clusters
library(ggraph)
clonalNetwork(seurat, 
              reduction = "umap", 
              identity = "ident",
              filter.identity = "2",
              cloneCall = "aa")

##高级桑吉图
library(circlize)
library(scales)
circles <- getCirclize(seurat, group.by = "new.cluster.ids")
#Just assigning the normal colors to each cluster
grid.cols <- scales::hue_pal()(length(unique(new.cluster.ids)))
names(grid.cols) <- levels(new.cluster.ids)
#Graphing the chord diagram
circlize::chordDiagram(circles, self.link = 1, grid.col = grid.cols)
ggsave("CD8 TCR sharing.PDF",width = 8,height = 8)


#不同患者和类群之间patient-cluster-type之间的关系，计算量大,反应的是转录组
seurat@meta.data$group[which(seurat$orig.ident == "P1"|seurat$orig.ident =="P5"|seurat$orig.ident =="P7"|seurat$orig.ident =="P9"|seurat$orig.ident =="P10")] = 'NR'
seurat@meta.data$group[which(seurat$orig.ident == "P2"|seurat$orig.ident =="P3"|seurat$orig.ident =="P4"|seurat$orig.ident =="P6"|seurat$orig.ident =="P8")] = 'R'
alluvialClonotypes(seurat, cloneCall = "gene", 
                   y.axes = c("orig.ident", "new.cluster.ids", "TRG"), 
                   color = "new.cluster.ids") 

#————————————————————————————————————————————————————————————————————
#提取扩增亚群并计算不同群之间的差异
df <- seurat[,seurat@meta.data$cloneType %in% "Hyperexpanded (100 < X <= 1000)" | seurat@meta.data$cloneType %in% "Large (20 < X <= 100)"]
#df <- seurat[,seurat@meta.data$cloneType %in% NA] #去除错误标记的细胞
df_R <- df[,df$orig.ident %in% "P3" | df$orig.ident %in% "P4" | df$orig.ident %in% "P6" | df$orig.ident %in% "P8"]
df_NR <- df[,df@meta.data$orig.ident %in% "P1" | df$orig.ident %in% "P10" | df$orig.ident %in% "P5" | df$orig.ident %in% "P7" | df$orig.ident %in% "P9"]
saveRDS(df,file="T cell_hyperexpanded & large.rds")
save(df_R,file = "T cell_R expand.rda")
save(df_NR,file = "T cell_NR expand.rda")
df <- get(load("T cell_hyperexpanded & large.rda"))
eso.freq = as.matrix(as.data.frame.matrix(table(df@meta.data$orig.ident,df@meta.data$CTaa)))
eso1.freq = as.matrix(as.data.frame.matrix(table(df@meta.data$orig.ident,df@meta.data$seurat_clusters)))
write.csv(eso.freq,"hyperexpanded.csv")
write.csv(eso1.freq,"expanded cluster.csv")

#提取相关信息进行STARTRAC--------------------------------------------------
meta = seurat@meta.data
NR = rbind(meta[which(meta$orig.ident == 'P1'),],meta[which(meta$orig.ident == 'P5'),],meta[which(meta$orig.ident == 'P7'),],meta[which(meta$orig.ident == 'P9'),],meta[which(meta$orig.ident == 'P10'),])
R = rbind(meta[which(meta$orig.ident == 'P2'),],meta[which(meta$orig.ident == 'P3'),],meta[which(meta$orig.ident == 'P4'),],meta[which(meta$orig.ident == 'P6'),],meta[which(meta$orig.ident == 'P8'),])
ids1 = NR[,c("new.cluster.ids","orig.ident","CTaa")] %>% as.data.frame() 
ids2 = R[,c("new.cluster.ids","orig.ident","CTaa")] %>% as.data.frame()
ids1$group = 'NR'
ids2$group = 'R'
df = rbind(ids1,ids2)
df = na.omit(df) #去除NA
df$loc = "T"
df=cbind(Cell_Name=row.names(df), df) #设置cell_name再去除rowname
row.names(df)=NULL
colnames(df)[2] = "majorCluster"
colnames(df)[3] = "patient"
colnames(df)[4] = 'clone.id'
write.table(df,"cloneDat_ESO2023_CD8.txt")

#三组比较
meta = seurat@meta.data
TRGI = rbind(meta[which(meta$orig.ident == 'P3'),],meta[which(meta$orig.ident == 'P4'),],meta[which(meta$orig.ident == 'P8'),])
TRGII = rbind(meta[which(meta$orig.ident == 'P2'),],meta[which(meta$orig.ident == 'P5'),],meta[which(meta$orig.ident == 'P6'),])
TRGIII = rbind(meta[which(meta$orig.ident == 'P1'),],meta[which(meta$orig.ident == 'P7'),],meta[which(meta$orig.ident == 'P9'),],meta[which(meta$orig.ident == 'P10'),])
ids1 = TRGI[,c("new.cluster.ids","orig.ident","CTaa")] %>% as.data.frame() 
ids2 = TRGII[,c("new.cluster.ids","orig.ident","CTaa")] %>% as.data.frame() 
ids3 = TRGIII[,c("new.cluster.ids","orig.ident","CTaa")] %>% as.data.frame() 
ids1$group = 'TRGI'
ids2$group = 'TRGII'
ids3$group ="TRGIII"
df = rbind(ids1,ids2,ids3)
df = na.omit(df) #去除NA
df$loc = "T"
df=cbind(Cell_Name=row.names(df), df) #设置cell_name再去除rowname
row.names(df)=NULL
colnames(df)[2] = "majorCluster"
colnames(df)[3] = "patient"
colnames(df)[4] = 'clone.id'
write.table(df,"cloneDat_ESO2023_CD8_TRG.txt")
