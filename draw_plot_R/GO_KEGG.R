#加载包
library(dplyr)
library(ggplot2)
library(readxl)
#输入go/kegg富集分析结果
a <- read_excel("~/Bio/后续富集分析流程/script/draw_plot_script/ZJst36F_vs_ZJst36M_go_all.xlsx")

#取前20个条目
b <- a %>%
  slice(1:20)

#GeneRatio 的数值转换
b$GeneRatioValue <- as.numeric(gsub("([^/]+)/.*", "\\1", b$GeneRatio)) /
  as.numeric(gsub(".*?/([^ ]+)", "\\1", b$GeneRatio))

# 绘制气泡图
ggplot(b, aes(x = GeneRatioValue, y = reorder(Description, Count), size = Count, color = pvalue)) +
  geom_point(shape = 16, alpha = 0.7) +  
  scale_color_gradient(low = "blue", high = "red", trans = "log10") +  # 颜色映射
  scale_size(range = c(2, 8)) +  
  labs(
    x = "Gene Ratio",
    y = "Description"
  ) +
  theme_bw() +  # 使用白色背景
  theme(
    axis.text.y = element_text(size = 10),
    legend.position = "right"             
  )
