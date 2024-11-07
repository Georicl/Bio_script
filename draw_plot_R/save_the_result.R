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
choose_org <- temp_args[3]
choose_exp <- temp_args[4]
#导入改名的sample
path_sample <- temp_args[5]
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

##################orgDB前期处理########################
#注意：该处理仅可以Ocu的包进行处理，同时会只提取第1、10、12、13、21列

#新库的构建方式，删除了重复行
emapper1 <- read_delim("~/Downloads/MM_j6g98g1i.emapper.annotations.tsv",    #读取原包orgDB
                       delim = "\t", escape_double = FALSE, col_names= FALSE,
                       comment = "#", trim_ws = TRUE)  %>%
  dplyr::select(GID=X1,GO=X10,KO=X12,Pathway=X13,OG=,Gene_name=X21) #删除一些不必要的列
emapper1$GID <-gsub("\\..*", "", as.character(emapper1$GID))  #删除末尾带.*的基因名
emapper1 <- emapper1 %>%  #删除重复行
  distinct(GID, .keep_all = TRUE) 
#删除GO列中包含'-'的行
empapper1 <- emapper1 %>%
  dplyr::filter(!grepl("-", GO))


write.table(empapper1,"~/Downloads/work/MM_xgaacda2.emapper.annotations.(1).txt",row.names = F,quote = F) #导出为处理后的库文件

####################orgDB处理后数据导入#################

#导入基因信息表，这个包是已经经过annotation处理的注释包
emapper <- read_table("~/Downloads/work/MM_xgaacda2.emapper.annotations.(1).txt")

#构建只含有GID和基因名的基因信息表
gene_info <- dplyr::select(emapper, GID,Gene_name) %>%
  dplyr::filter(!is.na(Gene_name))

#构建含有GID，GO，IEA的基因信息表
gene_GO <- dplyr::select(emapper, GID,GO)%>%
  mutate(EVIDENCE = 'IEA') %>%
  separate_rows(GO, sep = ',', convert = F)

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
save(gene_counts, gene_exp, genes_TMM_EXPR,
     sample_info,DE_result,
     file = './RNAseq.rdata')