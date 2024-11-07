library(cowplot)
library(ggplot2)
library(ggsci)
library(ggforce)
library(patchwork)
library(ggrepel)
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(clusterProfiler)
library(tidyverse)

load("./RNAseq.rdata")

############火山图绘图#################
#火山图
volcano_plot <- ggplot(DE_result, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(shape = 21, size = 4, aes(fill = ifelse(log2FoldChange < -1 & padj < 0.05, "down", 
                                                     ifelse(log2FoldChange > 1 & padj < 0.05, "up", "ns"))), 
             alpha = 0.8,) +  # 根据 log2FoldChange 值设置颜色，同时设定散点透明度与大小
  geom_text_repel(data = filter(DE_result, abs(log2FoldChange) > 1 & padj < 0.05),
                  box.padding = 0.3,
                  aes(label = id), #点名取自DE_result的id列
                  vjust = 1, hjust = 1, size = 3) +  # 控制标签位置和大小
  scale_fill_manual(
    values = c("down" = "blue", "up" = "red", "ns" = "gray")) +  #定义颜色
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +  # 添加 p 值参考线
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +  # 添加 FoldChange 阈值参考线
  theme_classic() +#使用经典绘图，没有背景图
  theme(
    panel.background = element_rect(fill = "transparent", color = NA),  # 绘图区域背景透明
    plot.background = element_rect(fill = "transparent", color = NA),   # 整个图形背景透明
  )
  labs(x = "Log2 Fold Change", y = "p-value",fill = 'significant')  #标题名，x轴名，y轴名
ggsave(file ="volcano_plot" ,width = 16, height = 12,bg = "transparent")
###########热图数据处理###################
#三表关联,用于制作热图
TNP <- dplyr::select(DE_result,id,log2FoldChange,pvalue,padj) %>% #选择出DE列，筛选出id，log2，pvalue，padj列
  mutate(direction=if_else(abs(log2FoldChange)<1 | padj>0.05,'ns',if_else(log2FoldChange >= 1 ,'up','down'))) %>% #设置不显著 上调和下调
  left_join(genes_TMM_EXPR, by = c('id'='id'))       #合并TMM表格和DE表格，用id列对照
#统计差异基因
group_by(TNP,direction) %>% 
  summarise(count=n())
#############总表达量热图##########
exp <-column_to_rownames(genes_TMM_EXPR,var = "id")
exp_sub <- exp[apply(exp, 1, function(x) sum(x) >= 1 ), ]
#指定文件名和尺寸
pdf(file = "Heatmap_expression_plot.png", width = 16, height = 12,bg = "transparent")
# 创建热图
Heatmap_expression_plot <- Heatmap(exp_sub,
                                   border = TRUE,
                                   name = "Expression",    # 设置颜色条名称
                                   col = colorRamp2(c(-5, 0, 5), c("blue", "white", "red")), # 设置渐变色
                                   show_row_names = TRUE,  # 显示行名
                                   show_column_names = TRUE) # 显示列名

# 绘制热图
draw(Heatmap_expression_plot)
dev.off()
## 计算相关性
dat_cor <- cor(t(exp_sub), method = "pearson")

#指定文件名和尺寸
pdf("Heatmap_Correlation_plot.png", width = 16, height = 12,bg = "transparent")
## 相关性热图
Heatmap_Correlation_plot <- pheatmap( dat_cor,
                                      color = colorRampPalette(c("blue","white","red"))(100),
                                      cluster_rows = F, 
                                      cluster_cols = F,
                                      display_numbers = F,
                                      # cellwidth = 10,
                                      # cellheight = 10
)
# 绘制热图
draw(Heatmap_Correlation_plot)
dev.off()
############组内差异热图#######
##热图
temp <- arrange(TNP,desc(abs(log2FoldChange))) %>%              #选择表格，以及用于制作热值的列           
  dplyr::select(-log2FoldChange,-pvalue, -padj,-direction) %>%    #选择计算的列
  column_to_rownames(var='id')  
temp <- log2(temp+0.001)
#temp转换为data.frame
temp_df <- as.matrix(temp)

#指定文件名和尺寸
pdf(file = "Heatmap_significant_plot.png", width = 16, height = 12,bg = "transparent")
#使用 ComplexHeatmap 绘制热图
Heatmap_significant_plot <-Heatmap(temp_df,
                                   name = "significant",    # 设置颜色条名称
                                   column_title = "Expression",
                                   col = colorRamp2(
                                     breaks = c(-5,0,5),
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
                                   column_km = 2,
                                   column_gap = unit(0.01,'npc'))
# 绘制热图
draw(Heatmap_significant_plot)  # 使用 draw() 函数来绘制热图
dev.off()  # 关闭图形设备，保存图像