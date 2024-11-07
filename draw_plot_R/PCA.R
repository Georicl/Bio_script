library(PCAtools)
library(ggalt)
library(readr)
library(BiocManager)
library(dplyr)
library(ggplot2)
library(tidyr)
library(tidyverse)
library(clusterProfiler)

###参数读取
sample_args <- commandArgs(trailingOnly = TRUE)
load('./RNAseq.rdata')
path_sample3 <- sample_args[1]
#########PCA
#样本信息表
sample2 <- read_delim(path_sample3, 
                      delim = "\t", escape_double = FALSE, 
                      col_names = FALSE, trim_ws = TRUE) 
sample_info2 <-rownames_to_column(sample_info)
#PCA前分组
gene_group <- sample2 %>%
  dplyr::select('X1','X2') %>%
  left_join(sample_info2, by = c('X2'='rowname'))
#给group组命名
gene_group1 <- gene_group %>%
  dplyr::rename(group =X1) %>%
  column_to_rownames(var = "group") %>%
  dplyr::rename(group = "X2")

#排序
gene_exp <- gene_exp[, order(colnames(gene_exp))]
gene_group1 <- gene_group1[order(rownames(gene_group1)), ]


#PCA主成分分析
PCA <- pca(gene_exp,
           metadata = gene_group1,
           removeVar = 0.1)
#载入pca数据
pca_loadings <- PCA[["loadings"]] #pca_loads表
pca_rotated  <- PCA[["rotated"]]  #载入rotated表
#可视化pca
##解释度,画图
screeplot(PCA)

##主成分差异,画图
biplot(PCA,x = "PC1",y = "PC2",
       encircle = TRUE,encircleFill = TRUE,
       colby = "group"
)

dev.off()  #关闭图形设备，保存图像