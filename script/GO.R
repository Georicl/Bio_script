library(org.Oc.eg.db,lib='~/Bio/orgDB/')
library(PCAtools)
library(ggalt)
library(readr)
#导入基因差异表达数据处理脚本中处理好的文本
###############数据导入###################
#准备样本信息表
sample_info <-read.table("~/Downloads/result/3.Merge_result/genes.TPM.not_cross_norm.TMM_info.txt",header = T, row.names = 1)

#准备表达矩阵，导入read counts数据
gene_counts <- read.table("~/Downloads/result/3.Merge_result/genes.counts.matrix",header = T, row.names = 1)

#导入表达矩阵
gene_exp <- read.table("~/Downloads/result/3.Merge_result/genes.TMM.EXPR.matrix",header = T, row.names = 1)

#导入表达矩阵，用于出图
genes_TMM_E <- read_delim("~/Downloads/result/3.Merge_result/genes.TMM.EXPR.matrix", 
                          delim = "\t", escape_double = FALSE, 
                          trim_ws = TRUE) %>%
  dplyr::rename(id=...1)
#DE结果
DE_r <- read_delim("~/Downloads/result/4.DE_analysis/mVSf/genes.counts.matrix.SYFGonad_vs_SYMGonad.DESeq2.DE_results", 
                   delim = "\t", escape_double = FALSE, 
                   trim_ws = TRUE)
DE_r <- DE_r %>%
  dplyr::rename(id=1)

#载入mapping对照表
sample3 <-read.table("~/Downloads/result/sample3.txt")
# 首先，将sample3的列名提取为命名向量
name_mapping <- setNames(sample3$V1, sample3$V2)

# 使用colnames()将gen_exp中的列名替换
colnames(gene_exp) <- ifelse(colnames(gene_exp) %in% names(name_mapping),
                             name_mapping[colnames(gene_exp)],
                             colnames(gene_exp))

# 使用colnames()将gen_counts中的列名替换
colnames(gene_counts) <- ifelse(colnames(gene_counts) %in% names(name_mapping),
                                name_mapping[colnames(gene_counts)],
                                colnames(gene_counts))
###########################数据处理#################

#导入基因信息表，这个包是已经经过annotation处理的注释包
emapper <- read_table("~/Downloads/work/MM_xgaacda2.emapper.annotations.(1).txt")

#构建只含有GID和基因名的基因信息表
gene_info <- dplyr::select(emapper, GID,Gene_name) %>%
  dplyr::filter(!is.na(Gene_name))

#筛选gene_exp的内容
columns_to_extract <- c(
  paste0(DE_r$sampleA, "1"),
  paste0(DE_r$sampleA, "2"),
  paste0(DE_r$sampleA, "3"),
  paste0(DE_r$sampleA, "4"),
  paste0(DE_r$sampleB, "1"),
  paste0(DE_r$sampleB, "2"),
  paste0(DE_r$sampleB, "3"),
  paste0(DE_r$sampleB, "4")
)
columns_to_extract <- unique(columns_to_extract)
# 从 gene_exp 提取对应列
gene_exp_f <- gene_exp %>%
  dplyr::select(any_of(columns_to_extract))

#DE结果处理
DE_r <- 
  # 添加上下调信息
  mutate(DE_r, direction = if_else(
    padj > 0.05, 'ns', if_else(
      abs(log2FoldChange) < 1, 'ns', if_else(
        log2FoldChange >= 1, 'up', 'down')
    ))) %>%
  # 关联基因信息
  left_join(gene_info, by = c('id' = 'GID')) %>%
  # 关联表达量信息
  left_join(rownames_to_column(gene_exp_f, var = 'id'), by = 'id') %>%
  # 去除无用的列
  dplyr::select(-c(2:6)) %>%
  # 按 log2FoldChange 绝对值降序排列
  arrange(desc(abs(log2FoldChange)))

genomic <- read.table("~/Downloads/result/genomic对照表.txt",header = TRUE) #选择对照表格


# 提取第一列和第二列
genomic_ids <- genomic[[1]]  # 第一列
new_names <- genomic[[1]]     # 第二列

# 创建一个数据框以便进行重命名
rename_df <- data.frame(old_id = genomic_ids, new_id = new_names)

# 从DE_result中保留与P450匹配的行
DE_result <- DE_r %>%
  filter(id %in% genomic_ids) %>%
  # 进行重命名
  left_join(rename_df, by = c("id" = "old_id")) %>%
  mutate(id = new_id) %>%    # 将id列更新为新名字
  dplyr::select(-new_id)             # 删除不再需要的列

# 定义要删除的列名列表
cols_to_remove <- c("SYMGonad1", "ZJDDMBrain1")

# 动态删除存在的列
DE_result <- DE_result %>%
  dplyr::select(-dplyr::any_of(cols_to_remove))



##############GO富集###############


#enrichGO富集分析模块
gene <-dplyr::filter(DE_result,abs(log2FoldChange)>2 & padj <0.05) %>%   #制作想要分析的基因列表，从DE_result中筛选
  pull('id')

EGO <- enrichGO(gene = gene,              #输入基因
                OrgDb = org.Oc.eg.db,     #输入构建的库
                keyType = "GID",          #数据类型
                ont = 'ALL',             #富集类型
                qvalueCutoff = 0.05,      #切断值
                pvalueCutoff = 0.05,      #切断值
)
#EGO输出后是一个结果文档，需要转换为数据框
ego_df <-as.data.frame(EGO)
#以padj排序并取前十作图
ego_top10<-ego_df %>%
  dplyr::arrange(p.adjust) %>%
  head(10)
ego_top10$Description<- gsub(",", "\n", ego_top10$Description)
#可视化
##柱状图
ggplot(ego_top10, aes(x = reorder(Description, Count), y = Count, fill = p.adjust)) +  #设置x轴和y轴，以及填充对象
  geom_col()+   #创建柱状图
  coord_flip() +    #横向展示
  scale_fill_gradient(low = "red", high = "blue", name = "p.adjust") +   #柱子填充的颜色，从低到高的颜色，并且根据p.adjust数值变化
  labs(x = "gene term", y = "Gene Count", title = "GO Enrichment") +  # XY轴，以及名字
  theme_minimal()    #简约主题

##箱线图
ggplot(ego_top10,aes(x= reorder(Description, Count),y = GeneRatio,fill = p.adjust))+ #xy轴
  geom_boxplot(aes(fill = GeneRatio),outlier.shape = 21)+         #箱线图
  coord_flip()+          #横向
  geom_jitter(width = 0.5,size =2)+
  labs(x = "gene term", y = "Gene Count", title = "GO Enrichment") +#添加散点
  theme_minimal()+
  theme(axis.text.y = element_blank())# 使用简约主题


##散点图
ggplot(ego_top10, aes(x = FoldEnrichment, y = Description,size = Count,colour = p.adjust))+    #选择x轴和y轴
  scale_color_gradient(low = "blue", high = "red")+
  geom_point()+  #散点图
  theme_bw()+   #设置主题
  
  scale_size(range = c(1,5)) #散点大小
#关于ggplot的可选项
# geom_bar(stat = "identity") 表示柱子自动识别提供的y值进行绘图，stat参数可以选择“count=统计每个类型的数量”
#enrich_result <- enrich_result %>%           #选择count的占比百分数
#mutate(Percent = Count / sum(Count) * 100)    #此时y =percent

#theme(axis.text.x = element_text(angle = 45, hjust = 1),  #旋转X轴标签
#plot.title = element_text(hjust = 0.5) )            # 标题居中
# theme_bw()：黑白主题，带有明显的网格线，适合展示科学数据。
#	theme_classic()：经典风格，只保留坐标轴和少量元素。
#	theme_void()：完全空白的主题，没有轴线和背景，适合自由设计。
#	theme_dark()：黑色背景，适合用于演示和数据视觉化。
# axis.title.x 和 axis.title.y：修改坐标轴标题样式。
#	axis.text.x 和 axis.text.y：修改坐标轴标签样式，例如字体大小、颜色、角度。
#	axis.line：修改坐标轴线的样式。
# legend.position：指定图例的位置，例如 "top", "bottom", "left", "right", "none".
#	legend.title：修改图例标题样式。（前面要加legend.title=element_text,说明自己要修改的是自定义文本的内容 如：face="bold"加粗）
#	legend.text：修改图例内容的样式。
# element有以下类型：line、rect、text （线条、矩形、文本）
# scale_fill_gradient：自定义填充颜色








