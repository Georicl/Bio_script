######### 设置参数 #######
library(getopt)
spec <- matrix(
  c("sample_path", "s", 2, "character", "input your sample_info_txt, maybe name like <genes.TPM.not_cross_norm.TMM_info.txt>",
    "count_path", "c", 2, "character", "gene count matrix",
    "expression_path", "e", 2, "character", "gene expression standardize matrix，check the first col must have name",
    "de_result", "d", 2, "character", "DE_result matrix, check the first col must have name",
    "map", "m", 1, "character", "input your mapping txt <gene>(the first colunm) map <transcript>(the second)",
    "output", "o", 2, "character", "output dir"),
    byrow = TRUE, ncol = 5
)

######## 解析参数 ######
opt_par <- getopt(spec = spec)

if (is.null(opt_par$sample_path) || is.null(opt_par$count_path) || is.null(opt_par$expression_path) || 
    is.null(opt_par$de_result) || is.null(opt_par$output)) {
  cat("Usage:\n")
  cat("Rscript script_name.R -s <sample_info_txt> -c <count_matrix> -e <expression_matrix> -d <de_result_matrix> -o <output_dir> [-m <mapping_txt>]\n\n")
  cat("Options:\n")
  cat("-s, --sample_path\tinput your sample_info_txt, maybe name like <genes.TPM.not_cross_norm.TMM_info.txt>\n")
  cat("-c, --count_path\tgene count matrix\n")
  cat("-e, --expression_path\tgene expression standardize matrix, check the first col must have header\n")
  cat("-d, --de_result\tDE_result matrix, check the first col must have header\n")
  cat("-m, --map\tinput your mapping txt <gene>(the first column) map <transcript>(the second)\n")
  cat("-o, --output\toutput dir\n")
  quit(status = 1)
}

######## 加载库 ######
library(tidyverse)

####### 解析参数 #####
sample_path <- opt_par$sample_path
count_path <- opt_par$count_path
expression_path <- opt_par$expression_path
de_result <- opt_par$de_result
output_path <- opt_par$output

if( !is.null(opt_par$map)){
  mapping_file <- opt_par$map
}

output_path_heatmap <- file.path(output_path, "heatmap_plot.pdf")
output_path_volcano <- file.path(output_path, "volcano.pdf")
####### 载入文件 #####
#准备样本信息表
sample_info <-read.table(sample_path,header = T, row.names = 1)

#准备表达矩阵，导入read counts数据
gene_counts <- read.table(count_path,header = T, row.names = 1)

#导入表达矩阵
gene_exp <- read.table(expression_path,header = T)

#导入表达矩阵，用于出图
genes_TMM_E <- read_delim(expression_path ,
                          delim = "\t", escape_double = FALSE, 
                          trim_ws = TRUE)
colnames(genes_TMM_E)[1] <- "id"

#DE结果导入，并且将第一列改名为id
DE_result <- read_delim(de_result, 
                        delim = "\t", escape_double = FALSE, 
                        trim_ws = TRUE)
colnames(DE_result)[1] <- "id"

####### 处理map参数 ###
if (!is.null(opt_par$map)) {
  colnames(gene_exp)[1] <- "id"
  # 读取 map 文件 (基因名与 DE_result 的 id 映射关系)
  gene_map <- read.table(mapping_file, header = FALSE, col.names = c("gene_name", "id"))
  
  # 确保 DE_result 和 gene_exp 的 id 列都是字符型
  DE_result$id <- as.character(DE_result$id)
  gene_exp$id <- as.character(gene_exp$id)
  genes_TMM_E$id <- as.character(genes_TMM_E$id)
  
  # 根据 map.txt 中的 gene_name 筛选 DE_result、gene_exp 和 genes_TMM_E
  DE_result <- DE_result %>%
    filter(id %in% gene_map$id)  # 只保留在 map.txt 中存在的 id
  
  gene_exp <- gene_exp %>%
    filter(id %in% gene_map$id)  # 同样筛选 gene_exp
  
  genes_TMM_E <- genes_TMM_E %>%
    filter(id %in% gene_map$id)  # 同样筛选 genes_TMM_E
  
  # 使用 gene_map 替换 DE_result 中的 id 为 gene_name
  DE_result <- DE_result %>%
    left_join(gene_map, by = "id") %>%
    mutate(id = gene_name) %>%
    select(-gene_name)  # 移除 gene_name 列
  
  gene_exp <- gene_exp %>%
    left_join(gene_map, by = "id") %>%
    mutate(id = gene_name) %>%
    select(-gene_name)  # 移除 gene_name 列
  
  genes_TMM_E <- genes_TMM_E %>%
    left_join(gene_map, by = "id") %>%
    mutate(id = gene_name) %>%
    select(-gene_name)  # 移除 gene_name 列


####### 数据处理 #####
# log2FC信息改成可选项
# DE结果处理
DE_result <- 
  # 添加上下调信息
  mutate(DE_result, direction = if_else(
    padj > 0.05, 'ns', if_else(
      abs(log2FoldChange) < 1, 'ns', if_else(
        log2FoldChange >= 1, 'up', 'down')
    ))) %>%
  # 关联表达量信息
  left_join(gene_exp, by = 'id') %>%  
  # 去除无用的列
  dplyr::select(-c(2:6,8,9))  %>%
  # 按 log2FoldChange 绝对值降序排列
  arrange(desc(abs(log2FoldChange)))
}else{
####### 数据处理(无map参数) #####
# DE结果处理
DE_result <- 
  # 添加上下调信息
  mutate(DE_result, direction = if_else(
    padj > 0.05, 'ns', if_else(
      abs(log2FoldChange) < 1, 'ns', if_else(
        log2FoldChange >= 1, 'up', 'down')
    ))) %>%
  # 关联表达量信息
  left_join(rownames_to_column(gene_exp, var='id'), by = 'id') %>%  
  # 去除无用的列
  dplyr::select(-c(2:6,8,9))  %>%
  # 按 log2FoldChange 绝对值降序排列
  arrange(desc(abs(log2FoldChange)))
}
################# 绘图 #################
library(ComplexHeatmap)
library(circlize)

## 热图
TNP <- dplyr::select(DE_result,id,log2FoldChange,pvalue,padj) %>%    #选择出DE列，筛选出id，log2，pvalue，padj列
  mutate(direction=if_else(abs(log2FoldChange)<1 | padj>0.05,'ns',if_else(log2FoldChange >= 1 ,'up','down'))) %>% #设置不显著 上调和下调
  left_join(genes_TMM_E, by = c('id'='id'))       #合并TMM表格和DE表格，用id列对照

#temp <- arrange(TNP,desc(abs(log2FoldChange))) %>%              #选择表格，以及用于制作热值的列  
#  dplyr::select(-log2FoldChange,-pvalue, -padj,-direction) %>%    #选择计算的列
#  column_to_rownames(var='id')
temp <- arrange(TNP, desc(abs(log2FoldChange))) %>%  # 按绝对值排序
  slice_head(n = 20) %>%                  # 取前20行
  dplyr::select(-log2FoldChange, -pvalue, -padj, -direction) %>% 
  column_to_rownames(var = 'id')      

temp <- log2(temp+1)
temp_df <- as.matrix(temp)

pdf(file = output_path_heatmap, width = 8, height = 6,bg = "transparent") 
p1 <- Heatmap(temp_df,
        name = "significant",    # 设置颜色条名称
        col = colorRamp2(
          breaks = c(-3,0,3),
          colors = c('blue','white','red')
        ),
        show_row_names = TRUE,  # 显示行名
        show_column_names = TRUE , # 显示列名
        border = 'grey',     #边框灰色
        rect_gp = gpar(col = 'white',lwd =1),   #内边框白色，宽度1
        row_names_gp = gpar(fontsize = 8,fontface = 'italic'), #字体
        column_names_gp = gpar(fontsize = 8,fontface = 'EUC'), #字体
        column_title_gp = gpar(fontsize = 16,fontface = 'italic'), #标题字体
        #聚类分组
        column_km =1 ,
        column_gap = unit(0.01,'npc'),
        use_raster = FALSE ,
)
draw(p1)
dev.off()
## 火山图
library(ggplot2)
library(ggsci)
library(cowplot)
library(ggrepel)
ggplot(DE_result, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(shape = 21, size = 4, aes(fill = ifelse(log2FoldChange < -1 & padj < 0.05, "down",
                                                     ifelse(log2FoldChange > 1 & padj < 0.05, "up", "ns"))),
             alpha = 0.8,) +  # 根据 log2FoldChange 值设置颜色，同时设定散点透明度与大小
  geom_text_repel(data = filter(DE_result, abs(log2FoldChange) > 3 & padj < 0.05),
                  box.padding = 0.3,
                  aes(label = id), #点名取自DE_result的id列
                  vjust = 1, hjust = 1, size = 3) +  # 控制标签位置和大小
  scale_fill_manual(
    values = c("down" = "blue", "up" = "red", "ns" = "gray")) +  #定义颜色
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +  # 添加 p 值参考线
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +  # 添加 FoldChange 阈值参考线
  theme_classic() +  #使用经典绘图，没有背景图
  labs(title = "volcano", x = "Log2 Fold Change", y = "p-value",fill = 'significant') + #标题名，x轴名，y轴名
  theme(plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "transparent", color = NA),  # 绘图区域背景透明
        plot.background = element_rect(fill = "transparent", color = NA),   # 整个图形背景透明
        legend.background = element_rect(fill = "transparent", color = NA) )  # 使得标题居中
ggsave(file = output_path_volcano, width = 8, height = 6,bg = "transparent")