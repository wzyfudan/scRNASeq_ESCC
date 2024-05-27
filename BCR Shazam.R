
rm(list=ls())
setwd("E:/Lab data/Sunlab/bioinformatics/rawdata")
#BiocManager::install("GenomicAlignments")
#The Change-O commandline tools provide a set of utilities for automated processing of 
  #Ig repertoire data following germline segment assignment using a tools such as IMGT/HighV-QUEST.
library(airr)      #载入适应性免疫组库
library(dplyr)

#计算SHM distance距离，change-O defineclone的cut-off阈值
#原数据中CDR3-junction,D-segemnt连接V和J
library(shazam)    
data <- read.table("SHM/heavy_parse-select.tsv", header=T, sep="\t")
matches <- data$sequence_alignment == data$germline_alignment
sum(matches)
#is_null_x <- is.null(data$c_call) | data$c_call == ""
#sum(is_null_x)，使用heavy_parse-select过滤结果去掉了空行

#获取changeO所需要的聚类阈值
# Use nucleotide Hamming distance and normalize by junction length
dist_ham <- distToNearest(data, sequenceColumn="junction", 
                          vCallColumn="v_call", jCallColumn="j_call",
                          model="ham", normalize="len", nproc=1)
#可视化,手动寻找阈值
library(ggplot2)
p1 <- ggplot(subset(dist_ham, !is.na(dist_nearest)),
             aes(x=dist_nearest)) + 
  theme_bw() + 
  xlab("Hamming distance") + 
  ylab("Count") +
  scale_x_continuous(breaks=seq(0, 1, 0.1)) +
  geom_histogram(color="white", binwidth=0.02) +
  geom_vline(xintercept=0.12, color="firebrick", linetype=2)
plot(p1)
# Find threshold using density method，自动输出阈值
#方法有 method有density或者gmm，前者为官网推荐的机器选择方式，后者为混合模型
output <- findThreshold(dist_ham$dist_nearest, method="density")
threshold <- output@threshold
# Plot distance histogram, density estimate and optimum threshold
plot(output, title="Density Method")
print(output)


#计算SHM使用函数“observedMutations”-----------------------------
library(shazam) 
data <- list()
{data[[1]] <- read.table("SHM/P1/heavy_parse-select_clone-pass_germ-pass.tsv", header=T, sep="\t")
data[[2]] <- read.table("SHM/P2/heavy_parse-select_clone-pass_germ-pass.tsv", header=T, sep="\t")
data[[3]] <- read.table("SHM/P3/heavy_parse-select_clone-pass_germ-pass.tsv", header=T, sep="\t")
data[[4]]<- read.table("SHM/P4/heavy_parse-select_clone-pass_germ-pass.tsv", header=T, sep="\t")
data[[5]] <- read.table("SHM/P5/heavy_parse-select_clone-pass_germ-pass.tsv", header=T, sep="\t")
data[[6]] <- read.table("SHM/P6/heavy_parse-select_clone-pass_germ-pass.tsv", header=T, sep="\t")
data[[7]] <- read.table("SHM/P7/heavy_parse-select_clone-pass_germ-pass.tsv", header=T, sep="\t")
data[[8]] <- read.table("SHM/P8/heavy_parse-select_clone-pass_germ-pass.tsv", header=T, sep="\t")
data[[9]] <- read.table("SHM/P9/heavy_parse-select_clone-pass_germ-pass.tsv", header=T, sep="\t")
data[[10]] <- read.table("SHM/P10/heavy_parse-select_clone-pass_germ-pass.tsv", header=T, sep="\t")}

{data[[1]]$sample <- "P1"
data[[2]]$sample <- "P2"
data[[3]]$sample <- "P3"
data[[4]]$sample <- "P4"
data[[5]]$sample <- "P5"
data[[6]]$sample <- "P6"
data[[7]]$sample <- "P7"
data[[8]]$sample <- "P8"
data[[9]]$sample <- "P9"
data[[10]]$sample <- "P10"}

{data[[1]]$group <- "NR"
data[[2]]$group <- "R"
data[[3]]$group <- "R"
data[[4]]$group <- "R"
data[[5]]$group <- "NR"
data[[6]]$group <- "R"
data[[7]]$group <- "NR"
data[[8]]$group <- "R"
data[[9]]$group <- "NR"
data[[10]]$group <- "NR"}

#批量计算mutation频率
seurat <- readRDS("B_cell_harmony.Rds")
meta <- seurat@meta.data
meta$sequence_id <- rownames(meta)
meta$sequence_id <- gsub('_.*',"",meta$sequence_id)
db <- list()
for (i in 1:10) {
db[[i]] <- observedMutations(
    data[[i]],
    sequenceColumn = "sequence_alignment",
    germlineColumn = "germline_alignment_d_mask",
    regionDefinition = NULL,
    mutationDefinition = NULL,
    ambiguousMode = c("eitherOr", "and"),
    frequency = T,
    combine = T,
    nproc = 1,
    cloneColumn = "clone_id",
    juncLengthColumn = "junction_length")  
db[[i]] <- db[[i]] %>%
  mutate(SHM = case_when(
    mu_freq > 0.03 ~ "high-SHM",
    mu_freq < 0.00001 ~ "no-SHM",
    TRUE ~ "low-SHM"))
db[[i]]$sequence_id <- gsub('_.*',"",db[[i]]$sequence_id)
db[[i]] <- merge(db[[i]],meta,by = "sequence_id")
#db[[i]] <- table(db[[i]][,c("sample","SHM")]) %>% as.data.frame()  #c_call表示IGMH等
#db[[i]] <- table(db[[i]][,c("new.cluster.ids","SHM","group.y")]) %>% as.data.frame()
#db[[i]]$prop = 100*db[[i]]$Freq/sum(db[[i]]$Freq)
}
#sum <- sum(db[[1]]$prop)
#计算频率
merged_table <- dplyr::bind_rows(db)  #合并db list
p<- ggplot(merged_table, aes(x = TRG, y = prop, fill = SHM)) +
  geom_bar(stat = "identity", position = "dodge") + facet_wrap(~ SHM, ncol = 1) + geom_point()+ylab("% among cells with productive BCR in Each Patient")
p
df_pvals <- merged_table %>%
  group_by(new.cluster.ids) %>% #横坐标
  rstatix::wilcox_test(  #选择t_test或者wilcox或者kruskal_test()
    prop ~ SHM,  #纵坐标[重点更改]和分组
    p.adjust.method = "BH",  #还可以选择bonferroni
    ref.group = "no-SHM") %>% 
  rstatix::add_significance("p") %>% #转化为
  rstatix::add_xy_position(x = "new.cluster.ids", dodge = 0.9)#确定p值的空间位置，否则会直接以分组显示
p + geom_text(data = merged_table[merged_table$SHM == "high-SHM", ], aes(x = 1, y = 75, label = "****"), color = "black", size = 5)+ 
  geom_text(data = merged_table[merged_table$SHM == "high-SHM", ], aes(x = 3, y = 75, label = "*"), color = "black", size = 5)+ 
  geom_text(data = merged_table[merged_table$SHM == "high-SHM", ], aes(x = 6, y = 75, label = "*"), color = "black", size = 5)+ 
  geom_text(data = merged_table[merged_table$SHM == "high-SHM", ], aes(x = 8, y = 75, label = "*"), color = "black", size = 5)+
  geom_text(data = merged_table[merged_table$SHM == "high-SHM", ], aes(x = 9, y = 75, label = "*"), color = "black", size = 5)
ggsave('cluster_SHM.PDF', width = 9, height = 6)

#计算发生了SHM的clone_id数量前面db()无需选取额外元素
clone <- list()
for (i in 1:10) {
clone[[i]] <- table(db[[i]][,c("clone_id","SHM","sample")]) %>% as.data.frame()
clone[[i]] <- clone[[i]] %>%
  group_by(clone_id) %>%
  mutate(total_SHM = sum(Freq),
         prop = Freq / total_SHM)
}
for (i in 1:10) {
  clone[[i]] <- subset(clone[[i]], Freq > 1 & SHM == "high-SHM" & prop > 0.5)
}
merge <- dplyr::bind_rows(clone)
merge <- table(merge[,c("sample")]) %>% as.data.frame()
merged_table <- dplyr::bind_rows(db)  #合并db list
merge1 <- merged_table %>%
  group_by(sample) %>%
  summarize(Freq = n_distinct(clone_id))
merge$group <- "R" 
merge$group[merge$sample %in% c("P1","P5","P7","P9","P10")] = "NR"
#{merge$TRG <- "TRG I" 
#merge$TRG[merge$sample %in% c("P2","P5","P6")] = "TRG II"
#merge$TRG[merge$sample %in% c("P1","P7","P9","P10")] = "TRG III"}

merge <- inner_join(merge, merge1, by = "sample")
merge <- merge %>%
  mutate(prop = Freq.x / Freq.y)
library(ggpubr)
ggplot(merge, aes(x = group, y = prop, color = group)) +
  geom_boxplot(width = 0.6) +
  ylab("B Clonetype with High-SHM proportion > 0.5 & Freq >1") +stat_compare_means(size = 5)+
  geom_point()+theme(axis.text.x = element_text(face = "bold", size = 15, hjust = 1))+
  stat_summary(fun = "mean",
               geom = "text",
               aes(label = round(..y.., 2)),
               position = position_dodge(width = 0.75),
               vjust = 2.3,
               size = 6)
ggsave("Patient clonetype SHM.PDF",width = 6,height = 6)
