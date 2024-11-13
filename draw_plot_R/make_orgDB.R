#加载必须库
library(BiocManager)
library(dplyr)
library(ggplot2)
library(tidyr)
library(tidyverse)
library(clusterProfiler)
library(pkgbuild)
library(AnnotationForge)
###参数读取
arg <- commandArgs(trailingOnly = TRUE)
path_orgdb <- arg[1] #这个参数要是orgDB的路径的
##################orgDB前期处理########################
#注意：该处理仅可以对Oryzias Curvintus的包进行了针对性优化，因为部分列为无用列，因此针对性的进行修改
#最终只提取第1、10、12、13、21列，并且会删除带*的基因名，如果orgDB导入的tsv符合如下规则，则可以直接导入
#第一列为GID
#第10列为GO编号
#第12列为KO编号
#第13列为Pathway通路
#第21列为Gene_name编号
#同时对于可能存在多个转录本的基因，改名后去除了转录变体，因此可能会出现重复行，所以需要去除重复行
emapper1 <- read_delim(path_orgdb,    #读取原包orgDB
                       delim = "\t", escape_double = FALSE, col_names= FALSE,
                       comment = "#", trim_ws = TRUE)  %>%
  dplyr::select(GID=X1,GO=X10,KO=X12,Pathway=X13,OG=,Gene_name=X21) #删除一些不必要的列
emapper1$GID <-gsub("\\..*", "", as.character(emapper1$GID))  #删除末尾带.*的基因名
emapper1 <- emapper1 %>%  #删除重复行
  distinct(GID, .keep_all = TRUE) 
# 删除 GO 列中包含 '-' 的行
empapper1 <- emapper1 %>%
  dplyr::filter(!grepl("-", GO))

#输出到当前目录
write.table(empapper1,"./MM_xgaacda2.emapper.annotations.(1).txt",row.names = F,quote = F) #导出为处理后的库文件
####################orgDB处理后数据导入#################

#导入基因信息表，这个包是已经经过annotation处理的注释包
emapper <- empapper1

#构建只含有GID和基因名的基因信息表
gene_info <- dplyr::select(emapper, GID,Gene_name) %>%
  dplyr::filter(!is.na(Gene_name))

#构建含有GID，GO，IEA的基因信息表
gene_GO <- dplyr::select(emapper, GID,GO)%>%
  mutate(EVIDENCE = 'IEA') %>%
  separate_rows(GO, sep = ',', convert = F)

##################orgDB构建，安装##########
#构建orgDB
AnnotationForge::makeOrgPackage(gene_info = gene_info,                     #基因信息表
                                go=gene_GO,             #GO基因分类
                                maintainer='man <walkinglight@163.com>',   #作者
                                author='1',      #作者
                                outputDir="./",    #输出目录
                                tax_id=0000,               #不重要
                                genus='O',                #属名
                                species='c',               #种名
                                goTable="go",             #GO识别
                                version="1.0")            #版本
# 打包
pkgbuild::build('./org.Oc.eg.db', dest_path = "./")

# 创建一个文件夹
dir.create('./orgDB', recursive = T)

# 将包安装在该文件夹下
install.packages('./org.Oc.eg.db_1.0.tar.gz'
                 , repos = NULL, lib='./orgDB')
