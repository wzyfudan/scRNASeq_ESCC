# source code downloaded from:
# https://github.com/cssmillie/ulcerative_colitis
# there is a 'ulcerative_colitis-master' folder which contains a number of r code scripts
rm(list=ls())
setwd("E:/Lab data/Sunlab/bioinformatics/rawdata")
## load required libraries (analysis.r code includes all the libraries for running all analyses,
## but as here only shows the Cellular Composition difference test, only libraries listed below need to be loaded.
library(Seurat)
library(RColorBrewer) #for brewer.pal
library(Matrix) #for Matrix
library(DirichletReg)
library(Formula)
library(data.table)
library(tidyverse)
library(cowplot)


## not all **.r in analysis.r need to be 'source'.
## only these three below are required for cellular composition difference test.
source('E:/Lab data/Sunlab/bioinformatics/bin/DR_test/mtx.r') #for sparse_cbind
source('E:/Lab data/Sunlab/bioinformatics/bin/DR_test/plot.r') #for matrix_barplot
source('E:/Lab data/Sunlab/bioinformatics/bin/DR_test/colors.r') #for set.colors
## this function is extracted from analysis.r 
dirichlet_regression = function(counts, covariates, formula){  
  # Dirichlet multinomial regression to detect changes in cell frequencies
  # formula is not quoted, example: counts ~ condition
  # counts is a [samples x cell types] matrix
  # covariates holds additional data to use in the regression
  #
  # Example:
  # counts = do.call(cbind, tapply(seur@data.info$orig.ident, seur@ident, table))
  # covariates = data.frame(condition=gsub('[12].*', '', rownames(counts)))
  # res = dirichlet_regression(counts, covariates, counts ~ condition)
  #ep.pvals = dirichlet_regression(counts=ep.freq, covariates=ep.cov, formula=counts ~ condition)$pvals
  
  # Calculate regression
  counts = as.data.frame(counts)
  counts$counts = DR_data(counts)
  data = cbind(counts, covariates)
  fit = DirichReg(counts ~ condition, data)
  
  # Get p-values
  u = summary(fit)
  #compared with healthy condition, 15 vars.
  pvals = u$coef.mat[grep('Intercept', rownames(u$coef.mat), invert=T), 4] 
  v = names(pvals)
  pvals = matrix(pvals, ncol=length(u$varnames))
  rownames(pvals) = gsub('condition', '', v[1:nrow(pvals)])
  colnames(pvals) = u$varnames
  fit$pvals = pvals
  
  fit
}

## Load metadata for discovery and validation cohorts
## data downloaded from: 
meta = read.csv('E:/Lab data/Sunlab/data/ESO scTCR-BCR/ESO sample information.txt', sep = "\t", header=T,fileEncoding="UTF-8") #csv-txt-UTF-8
# load seurat objects
eso.seur = readRDS("E:/Lab data/Sunlab/bioinformatics/rawdata/CD8T_cell.Rds")
eso.seur = readRDS("E:/Lab data/Sunlab/bioinformatics/rawdata/Mono_cell.Rds")
# set counts matrices
eso.counts = eso.seur@assays[['RNA']]@counts

# For the dirichlet-multinomial regression, we need to know the disease state for each sample
# We can get this from the metadata table as follows:
sample2health = data.frame(unique(data.frame(sample=gsub('', '', meta[,'sample']), health=meta[,'TRG'])), row.names=1)
sample2health = data.frame(unique(data.frame(sample=gsub('', '', meta[,'sample']), health=meta[,'response'])), row.names=1)
eso.cov = data.frame(condition=factor(sample2health[rownames(eso.freq),1], levels=c("TRG III",'TRG II',"TRG I")), row.names=rownames(eso.freq)) #以level第一个为默认control
eso.cov = data.frame(condition=factor(sample2health[rownames(eso.freq),1], levels=c("NR","R")), row.names=rownames(eso.freq)) 

# Calculate significant changes using dirichlet multinomial regression
# This returns a matrix of p-values for each cell type / disease state
# 3 condition:counts,cov,calculate P values
eso.pvals = dirichlet_regression(counts=eso.freq, covariates=eso.cov, formula=counts ~ condition)$pvals

# Plot cell proportions with u+s
eso.pct = 100*eso.freq/rowSums(eso.freq)
#visualization
p1 = matrix_barplot(eso.pct, group_by=eso.cov$condition, pvals=eso.pvals, colors=set.colors)
p1
save_plot(p1, file='Monocytes_proportion_TRG.pdf', nrow=1, ncol=1.5)
