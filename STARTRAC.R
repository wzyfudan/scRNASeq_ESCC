#https://github.com/jianye0383/STARTRAC
rm(list=ls())
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
setwd("E:/Lab data/Sunlab/bioinformatics/rawdata")
#devtools::install_github("jokergoo/ComplexHeatmap")
#devtools::install_github("Japrin/sscVis")
#devtools::install_github("Japrin/STARTRAC")
#BiocManager::install("impute")
#devtools::install_local("E:/Lab data/Sunlab/bioinformatics/bin/STARTRAC") #本地安装！code页面download，zip
library("Startrac")
library("tictoc")
library("ggpubr")
library("ggplot2")
library("ComplexHeatmap")
library("RColorBrewer")
library("circlize")

#数据读取The input to Startrac which the columns "cell_name",'clone.id', 'patient', 'majorCluster' and 'loc' are required.
in.dat <- read.table("article-PDF/cloneDat_ESO2023_CD8.txt",stringsAsFactors = F,head=T)
in.dat <- read.table("article-PDF/cloneDat_ESO2023_B.txt",stringsAsFactors = F,head=T)
#运行STARTRAC管道：
tic("Startrac.run")
out <- Startrac.run(in.dat[in.dat$group %in% c("R"),], proj="ESOtest", cores=NULL,verbose=F)
out1 <- Startrac.run(in.dat[in.dat$group %in% c("NR"),], proj="ESOtest", cores=NULL,verbose=F)
#导出表格数据，进行统计
out@cluster.data$group <- "R"
out1@cluster.data$group <- "NR"
df <- rbind(out@cluster.data,out1@cluster.data)
df <- df %>% 
  filter(aid != 'ESOtest') #去除ESOtest方便后续分析
#write.csv(df,"CD8_STARTRAC.csv")

#进行统计检验——expa作图
# t检验
#批量计算p值
library(rstatix)
#df$majorCluster <- factor(df$majorCluster, c("CXCL13 Tex", "IL7R Tm", "TNFSF9 Tcm",  "GZMK Tem", "CD74 Tem", "ITGAE Trm", "FCGR3A NK-T", "TRDV1 γδ-T", "CX3CR1 Temra","ISG15 Tex"))
#df$majorCluster <- factor(df$majorCluster, c("IL7R Tm", "LAIR2 Treg", "IFNG Th1/Tfh",  "RTKN2 Treg", "CXCL13 Th1/Tfh", "KLF2 Tm", "ISG15 Treg", "GZMA Tem", "TNFRSF9 Treg", "TNFRSF4 Treg"))

#比较cluster之间差异
df_pvals <- df %>%
  group_by(majorCluster) %>% #横坐标
  rstatix::wilcox_test(  #选择t_test或者wilcox
    expa ~ group,  #纵坐标和分组[[tran,expa]]
    p.adjust.method = "BH",  #还可以选择bonferroni
    ref.group = "NR") %>% 
  rstatix::add_x_position(x = "majorCluster", dodge = 0.9) #确定p值的空间位置
#改动x,y的标签
library(ggprism)
p <- NULL
p <- ggplot(df, aes(x=majorCluster,y=expa))+
  geom_boxplot(aes(fill = group))+
  theme_bw()
p <- p + add_pvalue(
  df_pvals, y = 0.4, xmin = "xmin", xmax = "xmax", tip.length = 0, 
  fontface = "italic", lineend = "round", bracket.size = 0.5
)
p
ggsave('B_expa_STARTRAC.PDF', width = 9, height = 6)
#比较grop之间
ggplot(df, aes(x=group,y=expa))+
  geom_boxplot(aes(fill = group))+
  theme_bw()+stat_compare_means()

ggsave('CD4_expa_STARTRAC_group.PDF', width = 6, height = 4)
ggsave('B_expa_NR_compare.PDF', width = 12, height = 9)

#trans作图
library(rstatix)
#df$majorCluster <- factor(df$majorCluster, c("CXCL13 Tex", "IL7R Tm", "TNFSF9 Tcm",  "GZMK Tem", "CD74 Tem", "ITGAE Trm", "FCGR3A NK-T", "TRDV1 γδ-T", "CX3CR1 Temra","ISG15 Tex"))
df$majorCluster <- factor(df$majorCluster, c("IL7R Tm", "LAIR2 Treg", "IFNG Th1/Tfh",  "RTKN2 Treg", "CXCL13 Th1/Tfh", "KLF2 Tm", "ISG15 Treg", "GZMA Tem", "TNFRSF9 Treg", "TNFRSF4 Treg"))
df_pvals <- df %>%
  group_by(majorCluster) %>% #横坐标
  rstatix::wilcox_test(  #选择t_test或者wilcox
    tran ~ group,  #纵坐标和分组
    p.adjust.method = "BH",  #还可以选择bonferroni
    ref.group = "NR") %>% 
  rstatix::add_x_position(x = "majorCluster", dodge = 0.9) #确定p值的空间位置
#改动x,y的标签
#点图
library(ggprism)
p <- NULL
p <- ggplot(df, aes(x=majorCluster,y=tran))+
  geom_boxplot(aes(fill = group))+
  theme_bw()
#增加P值,y调整P值位置
p <- p + add_pvalue(
  df_pvals, y = 1.6, xmin = "xmin", xmax = "xmax", tip.length = 0, 
  fontface = "italic", lineend = "round", bracket.size = 0.5
)
p
ggsave('B_trans_STARTRAC score.PDF', width = 12, height = 9)

#Gini-Simpson
library(rstatix)
#df$majorCluster <- factor(df$majorCluster, c("CXCL13 Tex", "IL7R Tm", "TNFSF9 Tcm",  "GZMK Tem", "CD74 Tem", "ITGAE Trm", "FCGR3A NK-T", "TRDV1 γδ-T", "CX3CR1 Temra","ISG15 Tex"))
#df$majorCluster <- factor(df$majorCluster, c("IL7R Tm", "LAIR2 Treg", "IFNG Th1/Tfh",  "RTKN2 Treg", "CXCL13 Th1/Tfh", "KLF2 Tm", "ISG15 Treg", "GZMA Tem", "TNFRSF9 Treg", "TNFRSF4 Treg"))
df_pvals <- df %>%
  group_by(majorCluster) %>% #横坐标
  rstatix::wilcox_test(  #选择t_test或者wilcox
    gini ~ group,  #纵坐标和分组
    p.adjust.method = "BH",  #还可以选择bonferroni
    ref.group = "NR") %>% 
  rstatix::add_x_position(x = "majorCluster", dodge = 0.9) #确定p值的空间位置
#改动x,y的标签
#点图
library(ggprism)
p <- NULL
p <- ggplot(df, aes(x=majorCluster,y=gini))+
  geom_boxplot(aes(fill = group))+
  theme_bw()
#增加P值,y调整P值位置
p <- p + add_pvalue(
  df_pvals, y = 1.6, xmin = "xmin", xmax = "xmax", tip.length = 0, 
  fontface = "italic", lineend = "round", bracket.size = 0.5
)
p
#比较grop之间
ggplot(df, aes(x=group,y=gini))+
  geom_boxplot(aes(fill = group))+
  theme_bw()+stat_compare_means()

ggsave('B_trans_STARTRAC score.PDF', width = 12, height = 9)





#三组比较——————————————————————————
in.dat <- read.table("cloneDat_ESO2023_CD8_TRG.txt",stringsAsFactors = F,head=T)
out1 <- Startrac.run(in.dat[in.dat$group %in% c("TRGI"),], proj="ESOtest", cores=NULL,verbose=F)
out2 <- Startrac.run(in.dat[in.dat$group %in% c("TRGII"),], proj="ESOtest", cores=NULL,verbose=F)
out3 <- Startrac.run(in.dat[in.dat$group %in% c("TRGIII"),], proj="ESOtest", cores=NULL,verbose=F)
#导出表格数据，进行统计
out1@cluster.data$group <- "TRGI"
out2@cluster.data$group <- "TRGII"
out3@cluster.data$group <- "TRGIII"
df <- rbind(out1@cluster.data,out2@cluster.data,out3@cluster.data)
df <- df %>% 
  filter(aid != 'ESOtest') #去除ESOtest方便后续分析

#比较p值
#df$majorCluster <- factor(df$majorCluster, c("CXCL13 Tex", "IL7R Tm", "TNFSF9 Tcm",  "GZMK Tem", "CD74 Tem", "ITGAE Trm", "FCGR3A NK-T", "TRDV1 γδ-T", "CX3CR1 Temra","ISG15 Tex"))
#df$majorCluster <- factor(df$majorCluster, c("IL7R Tm", "LAIR2 Treg", "IFNG Th1/Tfh",  "RTKN2 Treg", "CXCL13 Th1/Tfh", "KLF2 Tm", "ISG15 Treg", "GZMA Tem", "TNFRSF9 Treg", "TNFRSF4 Treg"))
df_pvals <- df %>%
  group_by(majorCluster) %>% #横坐标
  rstatix::wilcox_test(  #选择t_test或者wilcox或者kruskal_test()或fisher.test
    expa ~ group,  #纵坐标[重点更改]和分组
    p.adjust.method = "BH",  #还可以选择bonferroni
    ref.group = "TRGIII") %>%    #是否单侧less,greater
  rstatix::add_significance("p") %>% #转化为显著性
  rstatix::add_xy_position(x ="majorCluster", dodge = 0.8)  ####需要保留否则会错位,或者无法显示p值

#改动x,y的标签
#点图
library(ggprism)
library(ggpubr)
#增加P值
df_pvals <- subset(df_pvals, p < 0.2)
p <- NULL
p <- ggplot(df, aes(x=majorCluster,y=expa))+
  geom_boxplot(aes(fill = group))+
  theme(axis.text.x = element_text(face = "bold", 
                                     size = 12, 
                                     angle = 45, 
                                     hjust = 1))
p + stat_pvalue_manual(df_pvals,label = "p",tip.length = 0,y.position = c(0.75),size = 5)
ggsave('CD8T_expa_STARTRAC_TRG.PDF', width = 8, height = 6)


#比较grop之间
ggplot(df, aes(x=group,y=expa))+
  geom_boxplot(aes(fill = group))+
  theme_bw()+stat_compare_means(comparisons = list(c("TRGI","TRGIII"),c("TRGII","TRGIII")),size = 6) + theme(axis.text.x = element_text(size = 16))
ggsave('B_expa_STARTRAC_group_TRG.PDF', width = 10, height = 6)
