#加载必须库
library(BiocManager)
library(dplyr)
library(tidyr)
library(tidyverse)
library(clusterProfiler)
library(pkgbuild)
library(AnnotationForge)

#################参数数值导入###############
temp_args <- commandArgs(trailingOnly = TRUE)
path_temp <- temp_args[1] #导入了3.merge文件夹位置
path_DE_result <- temp_args[2] #导入DE
#两个可选参数，默认为FALSE
choose_exp <- temp_args[3]
#导入改名的sample
path_sample <- temp_args[4]
#导入改名的orgDB
path_orgdb <- temp_args[5] 
choose_org <- temp_args[6]
#解析文件夹内文件
path_sample_info <- file.path(path_temp,"genes.TPM.not_cross_norm.TMM_info.txt")
path_gene_counts <- file.path(path_temp,"genes.counts.matrix")
path_gene_exp <- file.path(path_temp,"genes.TMM.EXPR.matrix")
path_genes_TMM_EXPR <- file.path(path_temp,"genes.TMM.EXPR.matrix")

#####################导入数据###############

#准备样本信息表
sample_info <-read.table(path_sample_info,header = T, row.names = 1)

#准备表达矩阵，导入read counts数据
gene_counts <- read.table(path_gene_counts,header = T, row.names = 1)

#导入表达矩阵
gene_exp <- read.table(path_gene_exp,header = T, row.names = 1)

#导入表达矩阵，用于出图
genes_TMM_EXPR <- read_delim(path_genes_TMM_EXPR, 
                             delim = "\t", escape_double = FALSE, 
                             trim_ws = TRUE) %>%
  dplyr::rename(id=1)

#DE结果导入，并且将第一列改名为id
DE_result <- read_delim(path_DE_result, 
                        delim = "\t", escape_double = FALSE, 
                        trim_ws = TRUE)
DE_result <- DE_result %>%
  dplyr::rename(id=1)

##################if改名处理#########################
#区分是否改动了exp
if(choose_exp == 1){
    #读取sample文件
    sample3 <- read.table(path_sample,header = FALSE)
    #将sample3的列名提取
    name_mapping <- setNames(sample3$V1, sample3$V2)

    #将genes_TMM_EXPR中的列名替换
    colnames(genes_TMM_EXPR) <- ifelse(colnames(genes_TMM_EXPR) %in% names(name_mapping),
                                   name_mapping[colnames(genes_TMM_EXPR)],
                                   colnames(genes_TMM_EXPR))

    #将DE_results的列替换
    colnames(DE_result) <- ifelse(colnames(DE_result) %in% names(name_mapping),
                              name_mapping[colnames(DE_result)],
                              colnames(DE_result))
    #将gene_exp的列替换
    colnames(gene_exp) <- ifelse(colnames(gene_exp) %in% names(name_mapping),
                              name_mapping[colnames(gene_exp)],
                              colnames(gene_exp))
    #使用colnames()将gen_counts中的列名替换
    colnames(gene_counts) <- ifelse(colnames(gene_counts) %in% names(name_mapping),
                             name_mapping[colnames(gene_counts)],
                             colnames(gene_counts))

}

#####################数据处理######################

#DE结果处理
DE_result <- DE_result %>%
  #添加上下调信息
  mutate( direction = if_else(
    padj > 0.05, 'ns', if_else(
      abs(log2FoldChange) < 1, 'ns', if_else(
        log2FoldChange >= 1, 'up', 'down')
    ))) %>%
  #关联表达量信息
  left_join(rownames_to_column(gene_exp, var = 'id'), by = 'id') %>%
  #去除无用的列
  dplyr::select(-c(2:6,8,9)) %>%
  #按 log2FoldChange 绝对值降序排列
  arrange(desc(abs(log2FoldChange)))
#区分是否添加orgDB
if(choose_org == 1){
DE_result <- DE_result %>%
  left_join(gene_info, by = c('id' = 'GID')) 
}

#########保存数据#########################
#保存已经导入的数据
#这个数据会在运行了画图plot后删除，无需变更，只是一个临时的中间文件
save(gene_counts, gene_exp, genes_TMM_EXPR,
     sample_info,DE_result,
     file = './RNAseq.rdata')