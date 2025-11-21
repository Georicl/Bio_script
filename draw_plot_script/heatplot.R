##########################
# 该脚本接受两个输入参数，用于输出结果图片
# 参数1为表达矩阵，RNA-seq标准流程生成的名字如genes.TMM.EXPR.matrix
# 参数2为DE_result/egdeR的结果count
# 参数3设置为log2的阈值
# 参数4设置为p-value的阈值
# 参数5设置为heatmap量谱最小值
# 参数6设置为heatmap量谱最大值
# 需要两个导入表格中，第一列有对应名字，或者一个空格
# 作者：XiangYang
##########################

library(tidyverse)
library(readxl)
args <- commandArgs()

exp_table <- args[1]
de_table <- args[2]
log2fold <- args[3]
p_value <- args[4]
heatmap_min <- args[5]
heatmap_max <- args[6]
#导入表达矩阵
gene_exp <- read.table(exp_table, header = T)


#DE结果导入，并且将第一列改名为id
DE_result <- read.table(de_table, header = T)
colnames(DE_result)[1] <- "id"

################# 绘图 #################
library(ComplexHeatmap)
library(circlize)

# 热图
temp <- log2(gene_exp+1)
temp_df <- as.matrix(temp)

Heatmap(temp_df,
              name = "significant",    # 设置颜色条名称
              col = colorRamp2(
                breaks = c(heatmap_min, 0, heatmap_max),
                colors = c('blue', 'white', 'red')
              ),
              show_row_names = TRUE,  # 显示行名
              show_column_names = TRUE , # 显示列名
              border = 'grey',     #边框灰色
              rect_gp = gpar(col = 'white', lwd =1),   #内边框白色，宽度1
              row_names_gp = gpar(fontsize = 8, fontface = 'italic'), #字体
              column_names_gp = gpar(fontsize = 8, fontface = 'EUC'), #字体
              column_title_gp = gpar(fontsize = 16, fontface = 'italic'), #标题字体
              #聚类分组
              column_km =1 ,
              column_gap = unit(0.01,'npc'),
              use_raster = FALSE ,
)
## 火山图
library(ggplot2)
ggplot(DE_result, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(shape = 21, size = 4, aes(fill = ifelse(log2FoldChange < log2fold & padj < p_value, "down",
                                                     ifelse(log2FoldChange > log2fold & padj < p_value, "up", "ns"))),
             alpha = 0.8,) +  # 根据 log2FoldChange 值设置颜色，同时设定散点透明度与大小
  geom_text_repel(data = filter(DE_result, abs(log2FoldChange) > log2fold & padj < p_value),
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
