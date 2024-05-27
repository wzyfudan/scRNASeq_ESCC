rm(list=ls())
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(ggpubr)
setwd("E:/Lab data/Sunlab/bioinformatics/rawdata")
ESOtest <- readRDS("CD8T_cell.Rds")

#group
meta = ESOtest@meta.data
NR = rbind(meta[which(meta$orig.ident == 'P1'),],meta[which(meta$orig.ident == 'P5'),],meta[which(meta$orig.ident == 'P7'),],meta[which(meta$orig.ident == 'P9'),],meta[which(meta$orig.ident == 'P10'),])
R = rbind(meta[which(meta$orig.ident == 'P2'),],meta[which(meta$orig.ident == 'P3'),],meta[which(meta$orig.ident == 'P4'),],meta[which(meta$orig.ident == 'P6'),],meta[which(meta$orig.ident == 'P8'),])
ids1 = table(NR[,c("new.cluster.ids","orig.ident")]) %>% as.data.frame() 
ids2 = table(R[,c("new.cluster.ids","orig.ident")]) %>% as.data.frame()
ids1$group = 'NR'
ids2$group = 'R'
df = rbind(ids1,ids2)
df$Freq[df$Freq == 0]<-NA
df = na.omit(df)
# prop analysis
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
  group_by(new.cluster.ids) %>% 
  rstatix::wilcox_test(  
    prop ~ group,  
    p.adjust.method = "BH", 
    ref.group = "NR") %>% 
  rstatix::add_significance("p") %>% 
  rstatix::add_x_position(x = "new.cluster.ids", dodge = 0.9)
#write.csv(df_pvals,"df_pvals.csv")
#df_pvals <- read.csv("df_pvals.csv",row.names = 1)

#
library(ggprism)
library(ggpubr)

#df_pvals <- subset(df_pvals, p < 0.05)

library(ggplot2)
ggplot(df, aes(x=new.cluster.ids,y=prop))+geom_boxplot(aes(fill = group))+ 
  stat_pvalue_manual(df_pvals,label = "p",tip.length = 0.01,y.position = 0.3,size = 5)+ theme(axis.text.x = element_text(face = "bold", 
                                                                                                                          size = 12, 
                                                                                                                          angle = 45, 
                                                                                                                          hjust = 1))
ggsave("CD8T_wilcox_R.PDF",width = 8,height = 6)

#TRG analysis------------------------------------------------------------------------------------------------------------
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

#prop
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
  group_by(new.cluster.ids) %>% 
  rstatix::wilcox_test(  
    prop ~ group,  
    p.adjust.method = "BH",  
    ref.group = "TRGIII") %>%   
  rstatix::add_significance("p") %>% 
  rstatix::add_x_position(x = "new.cluster.ids", dodge = 0.8)  

#write.csv(df_pvals,"Mono_df_pvals.csv")
#df_pvals <- read.csv("Mono_df_pvals.csv",row.names = 1)

library(ggprism)
library(ggpubr)

df_pvals <- subset(df_pvals, p < 0.3)
p + stat_pvalue_manual(df_pvals,label = "p",tip.length = 0.01,y.position = c(0.2,0.2,0.2,0.18),size = 5) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))+ theme(axis.text.x = element_text(face = "bold", 
                                                                                                size = 12, 
                                                                                                angle = 45, 
                                                                                                hjust = 1))
ggsave('Mono_wilcox_TRG_all.pdf', width = 8, height = 6)

