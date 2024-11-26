---
title: "Readme"
output: html_notebook
author: 向旸
---
# **seqtools**
It integrates some functions for reading blast and hisat2 result files for screening and summarizing results, and supports extracting sequences by geneid, sequencing sequences, and randomly extracting or disrupting sequence order for comparison analysis.\
Requirement environment: python>3.5\
The tool relies on the Bio module and pysam module need to download\
How to check out the help instructions:`python seqtools.py `\
Now seqtools has eight functional modules\
`sort` Sorting the length of nucleotide sequences，the default parameters are shortest to longest,but --reverse option will reverses the result.\
`shuffle` Break up the order of the sequence.\
`subseq` Screening and extraction of gene sequences based on provided gene ids,need to provide a sample list,the first list must same to geneID。
`extract`Randomly extracted sequences, seed number defaults to 10, with the option to enter your own seed (but it seems like I haven't gotten around to writing it yet)\
`N50` This option can read the fasta input and create a fai and output N50 and total length and GC content in a sample.txt\
`get_longest` can get the longest transcripts in cds sequence by gff3\

# **draw_plot**
This script (draw_plot fold) set is to serve some data generated after the upstream analysis of the server cannot be quickly automated in the downstream at one time, so a simple script to summarize the call is also written as a python practice. There may be some other updates that will continue to be based on this framework.\
Requirement environment: python>3.5\
Please watch check out the help instructions：`python /seqtools/draw_plot.py --help`\
How to make your sample.txt?\
You need to review to your genes.TMM.EXPR.matrix file.\
Second colunm in sample.txt must same to the genes.TMM.EXPR.matrix file's colunm (I think the file's colunm name should be same to your sample's name), The first colunm is the group's name.\
If you need to build a orgDB package, please download annotation package by tsv format.\
This script's default output path and temporary file storage path are the current directory.\
log named draw_plot.log, save in the current directory.\
Script Run Logic：\
parse arguements --confirm theavailability of orgDB\
if not input orgDB pacakage, skip make orgDB and run save_the_results.R input variance analysis files,output the R data temporary file\
Use this rdata to draw variance analysis plot like volcano plot(use ggplot2) ,heatmap(ComplexHeatmap),and PCA cluster analysis.\
Georicl@outlook.com\
Xiang Yang\


# **draw_plot**
这个脚本集是为了服务于一些进行了服务器上游分析后，产生的数据无法一次性在下游快速自动化处理，因此简单写了一个总结调用的脚本，也当是写python练手。后续可能会更新一些别的根据这个框架继续更新的。\
需求环境：python>3.5\
查看帮助文档：`python draw_plot.py --help`\
sample.txt需要先看你的genes.TMM.EXPR.matrix文件\
第二列的名字要和上述文件的列名相同（也就是样本的名称），第一列则是组别。绝大多数原始文件的TMM.EXPR文件的样本名也是原始组名。\
如果需要设置改名参数，请自己准备一个txt，第一列为组别，第二列则为样本名。\
orgDB参数导入时，请下载tsv格式的annotation包，否则无法进行。\
该脚本默认输出路径与临时文件存放路径为当前目录。\
运行日志名为draw_plot.log，存放在当前目录。\
Georicl@outlook.com\
向旸\

# **seqtools**
整合了一些用于阅读blast和hisat2的结果文件用于筛选和结果汇总的功能，并支持按geneid提取序列，为序列排序，随机提取或打乱序列顺序用于比对分析等功能。\
需求环境：python>3.5\
该工具依赖于Bio模块和pysam模块需要下载\
查看帮助文档：`python /seqtools/seqtools.py `\
目前seqtools有八个功能模块\
`sort` 对核苷酸序列的长度进行排序，默认参数为最短到最长，但--reverse选项将反转结果。\
`shuffle` 打乱序列顺序。\
`subseq` 根据提供的基因id进行基因序列筛选和提取，需要提供一个样本列表，第一个列表必须与基因ID相同。\
`extract` 随机提取序列，种子数默认为10，但可以通过选项输入自己的种子（但我还没来得及写）。\
`N50` 此选项可以读取fasta输入并创建fai文件，输出N50、总长度和GC含量。\
`get_longest` 可以根据gff3文件获取cds序列的最长转录本。\