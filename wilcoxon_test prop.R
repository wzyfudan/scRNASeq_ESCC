rm(list=ls())
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(ggpubr)
setwd("E:/Lab data/Sunlab/bioinformatics/rawdata")
ESOtest <- readRDS("CD8T_cell.Rds")
#ESOtest <- ESOtest[,ESOtest@meta.data$seurat_clusters %in% c(0,2,3,4,5,9)]
#ESOtest <- ESOtest[,ESOtest@meta.data$seurat_clusters %in% c(0,6,7,8,5,9)]
#ESOtest@meta.data$kind[ESOtest$seurat_clusters %in% c(2,3,4)] <- "GZMK T"
#ESOtest@meta.data$kind[ESOtest$seurat_clusters %in% c(6,7,8)] <- "GZMB T"
#ESOtest@meta.data$kind[ESOtest$seurat_clusters %in% c(0,5,9)] <- "Tex"
#细胞比例在不同分组间的差异https://blog.csdn.net/qq_40966210/article/details/120804976
# 分组并重命名dataframe
meta = ESOtest@meta.data
NR = rbind(meta[which(meta$orig.ident == 'P1'),],meta[which(meta$orig.ident == 'P5'),],meta[which(meta$orig.ident == 'P7'),],meta[which(meta$orig.ident == 'P9'),],meta[which(meta$orig.ident == 'P10'),])
R = rbind(meta[which(meta$orig.ident == 'P2'),],meta[which(meta$orig.ident == 'P3'),],meta[which(meta$orig.ident == 'P4'),],meta[which(meta$orig.ident == 'P6'),],meta[which(meta$orig.ident == 'P8'),])
ids1 = table(NR[,c("new.cluster.ids","orig.ident")]) %>% as.data.frame() 
ids2 = table(R[,c("new.cluster.ids","orig.ident")]) %>% as.data.frame()
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
  group_by(new.cluster.ids) %>% #横坐标
  rstatix::wilcox_test(  #选择t_test或者wilcox或者kruskal_test()
    prop ~ group,  #纵坐标[重点更改]和分组
    p.adjust.method = "BH",  #还可以选择bonferroni
    ref.group = "NR") %>% 
  rstatix::add_significance("p") %>% #转化为
  rstatix::add_x_position(x = "new.cluster.ids", dodge = 0.9)#[重点更改]确定p值的空间位置，否则会直接以分组显示或是p不显示
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
ggplot(df, aes(x=new.cluster.ids,y=prop))+geom_boxplot(aes(fill = group))+ 
  stat_pvalue_manual(df_pvals,label = "p",tip.length = 0.01,y.position = 0.3,size = 5)+ theme(axis.text.x = element_text(face = "bold", 
                                                                                                                          size = 12, 
                                                                                                                          angle = 45, 
                                                                                                                          hjust = 1))
ggsave("CD8T_wilcox_R.PDF",width = 8,height = 6)

#三分类组------------------------------------------------------------------------------------------------------------
ESOtest <- readRDS("B_cell_harmony.rds")
meta = ESOtest@meta.data
TRGI = rbind(meta[which(meta$orig.ident == 'P3'),],meta[which(meta$orig.ident == 'P4'),],meta[which(meta$orig.ident == 'P8'),])
TRGII = rbind(meta[which(meta$orig.ident == 'P2'),],meta[which(meta$orig.ident == 'P5'),],meta[which(meta$orig.ident == 'P6'),])
TRGIII = rbind(meta[which(meta$orig.ident == 'P1'),],meta[which(meta$orig.ident == 'P7'),],meta[which(meta$orig.ident == 'P9'),],meta[which(meta$orig.ident == 'P10'),])
ids1 = table(TRGI[,c("new.cluster.ids","orig.ident")]) %>% as.data.frame() 
ids2 = table(TRGII[,c("new.cluster.ids","orig.ident")]) %>% as.data.frame() 
ids3 = table(TRGIII[,c("new.cluster.ids","orig.ident")]) %>% as.data.frame() 
ids1$group = 'TRGI'
ids2$group = 'TRGII'
ids3$group ="TRGIII"
df = rbind(ids1,ids2,ids3)
df$Freq[df$Freq == 0]<-NA
df = na.omit(df)
#colnames(df)[3] <- "clonality"
#计算prop
ESOtest1 = readRDS("scRNA_harmony_UMAP.Rds")
meta = ESOtest1@meta.data
TRGI = rbind(meta[which(meta$orig.ident == 'P3'),],meta[which(meta$orig.ident == 'P4'),],meta[which(meta$orig.ident == 'P8'),])
TRGII = rbind(meta[which(meta$orig.ident == 'P2'),],meta[which(meta$orig.ident == 'P5'),],meta[which(meta$orig.ident == 'P6'),])
TRGIII = rbind(meta[which(meta$orig.ident == 'P1'),],meta[which(meta$orig.ident == 'P7'),],meta[which(meta$orig.ident == 'P9'),],meta[which(meta$orig.ident == 'P10'),])
ids1 = table(TRGI[,c("orig.ident")]) %>% as.data.frame() 
ids2 = table(TRGII[,c("orig.ident")]) %>% as.data.frame()
ids3 = table(TRGIII[,c("orig.ident")]) %>% as.data.frame()
df1 = rbind(ids1,ids2,ids3)
colnames(df1)[1] <- "orig.ident"
merged_df <- left_join(df, df1, by = "orig.ident")
df$prop = merged_df[,3]/merged_df[,5]
p <- ggplot(data=df, mapping=aes(x=new.cluster.ids,y=prop,color=group))+
  geom_boxplot()+theme_bw()
p
df_pvals <- df %>%
  group_by(new.cluster.ids) %>% #横坐标
  rstatix::wilcox_test(  #选择t_test或者wilcox或者kruskal_test()或fisher.test
    prop ~ group,  #纵坐标[重点更改]和分组
    p.adjust.method = "BH",  #还可以选择bonferroni
    ref.group = "TRGIII") %>%    #是否单侧less,greater
  rstatix::add_significance("p") %>% #转化为显著性
  rstatix::add_x_position(x = "new.cluster.ids", dodge = 0.8)  #需要保留否则会错位

#write.csv(df_pvals,"Mono_df_pvals.csv")
#df_pvals <- read.csv("Mono_df_pvals.csv",row.names = 1)

#改动x,y的标签
#点图
library(ggprism)
library(ggpubr)
#增加P值
df_pvals <- subset(df_pvals, p < 0.3)
p + stat_pvalue_manual(df_pvals,label = "p",tip.length = 0.01,y.position = c(0.2,0.2,0.2,0.18),size = 5) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))+ theme(axis.text.x = element_text(face = "bold", 
                                                                                                size = 12, 
                                                                                                angle = 45, 
                                                                                                hjust = 1))
ggsave('Mono_wilcox_TRG_all.pdf', width = 8, height = 6)



#再比值----------------------------------------------
result <- df %>%
  group_by(orig.ident) %>%
  summarise(ratio = sum(prop[kind == "GZMK T"]) / sum(prop[kind == "Tex"]))
{result$group <- "R" 
  result$group[result$orig.ident %in% c("P1","P5","P7","P9","P10")] = "NR"}
library(ggpubr)
ggplot(result, aes(x=group,y=ratio))+geom_boxplot(aes(fill = group))+ 
  theme_bw()+ stat_compare_means()+ylab("GZMB/Tex ratio")
ggsave("GZMB_R.pdf",width = 6,height = 6)
{result$TRG <- "TRG I" 
 result$TRG[result$orig.ident %in% c("P2","P5","P6")] = "TRG II"
 result$TRG[result$orig.ident %in% c("P1","P7","P9","P10")] = "TRG III"}
ggplot(result, aes(x=TRG,y=ratio))+geom_boxplot(aes(fill = TRG))+ 
  theme_bw()+ stat_compare_means(comparisons = list(c("TRG I","TRG III"),c("TRG II","TRG III")))+ylab("GZMB/Tex ratio")
ggsave("GZMB_TRG.pdf",width = 6,height = 6)
