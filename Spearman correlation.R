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

df1 = rbind(ids1,ids2) #换个细胞合并再来一次
df <- rbind(df,df1)
#去除值为0的行
df$Freq[df$Freq == 0]<-NA
df = na.omit(df)
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
#报错重启
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




#pearson相关性图
scRNA_epi <- ESOtest1[,ESOtest1@meta.data$seurat_clusters %in% c(10)]
scRNA_epi <- AddModuleScore(object = scRNA_epi, features = "S100A7", name = "cell")
data = scRNA_epi@meta.data[,c("orig.ident","cell1","group")] %>% as.data.frame()
data <- aggregate(data$cell1, by=list(data$orig.ident),mean)
colnames(data)[1] <- "orig.ident"

#或者结合bulk-seq结果
data <- data.frame(orig.ident = c("P1","P2","P5","P6","P7","P8","P9","P10"),value = as.numeric(c("66.19577617",	"2.101928394",	"2.959163305",	"1.03356623",'308.8712722',"159.3664925","2528.581834","1198.500916")))
data <- data.frame(orig.ident = c("P1","P2","P5","P6","P7","P8","P9","P10"),value = as.numeric(c("3.597020147","0","0","0.034532703","2.72861569",	"0.076982575","40.01423987","19.28345611")))

data <- read.csv("E:/lab data/Sunlab/data/20230510-mRNA-WZY/5.Expression/All.tpm.exp_anno.csv",row.names = 1)
{TLS <- c("CCL2", "CCL3", "CCL4", "CCL5", "CCL8", "CCL18", "CCL19", "CCL21", "CXCL9","CXCL10", "CXCL11", "CXCL13") #在scRNA打分需要加上list
sub_data <- data[data$GeneSymbol %in% TLS,]
sub_data <- sub_data[,1:8] %>% data.frame() #手动计算上述基因集平均值
data <- data.frame(orig.ident = c("P1","P2","P5","P6","P7","P8","P9","P10"),value = as.numeric(c("35.15950717",	"35.87923833","5.225784692","311.2459892","38.95326333","29.80269368","27.674503","41.5071495")))}

data <- mean_data #根据addmodulescore得到，因为肿瘤细胞无关TLS形成
colnames(data)[2] <-"value"
#data <- data.frame(orig.ident = c("P1","P2","P3","P4","P5","P6","P7","P8","P9","P10"),value = as.numeric(c("0.32",	"-0.07","0.44","0.34","0.38","0.52","-0.04","1.28","0.27","0.37"))) #FORL2 Mean


sub_df <- df[df$new.cluster.ids %in% c("RTKN2 Treg","ISG15 Treg","LAIR2 Treg","TNFRSF4 Treg","TNFRSF9 Treg"),]  #结合某细胞的比例
result <- sub_df %>% #计算多个cluster的总占比
  group_by(orig.ident) %>%
  mutate(all = sum(prop))

#合并计算相关性
result <- merge(result,data,by = "orig.ident")
sub_df$value[sub_df$value == 0]<-NA
sub_df = na.omit(sub_df)
correlation <- cor.test(result$all,result$value, method = "pearson")
ggplot(data = result, aes(x = value, y = all)) +
  geom_point()
p <- ggplot(data = result, aes(x = all, y = value)) +
  geom_point() +
  labs(title = paste("Pearson Correlation =", round(correlation$estimate, 2), " (P =", format(correlation$p.value, digits = 2), ")"),
       x = "IGHG PC & IGHA PC",
       y = "12-chemokines for Tertiary lymphoid structure(TLS) Mean TPM") 
p +geom_smooth(method = "lm", se = TRUE, color = "red") 
ggsave("TLS & IGHG PC.PDF",width = 6,height = 5)
