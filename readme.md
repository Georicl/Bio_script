# Bio Script
Yang \
These scripts were written when I was conducting research.

# **CQtools**
该整合重测序比对后进行染色体商（CQ）计算流程，包含cqmapping比对流程和cqtools主计算流程\
主要依赖软件：\
samtools\
bedtools\
bwa\
更多帮助请寻求：`python main.py cqtools -h`\

# **seqtools**
It integrates some functions for reading blast and hisat2 result files for screening and summarizing results, and supports extracting sequences by geneid, sequencing sequences, and randomly extracting or disrupting sequence order for comparison analysis.\
Requirement environment: python>3.5\
The tool relies on the Bio module and pysam module need to download\
For linux software dependencies: Blast、samtools. \
How to check out the help instructions:`python main.py `\
Now seqtools has seven functional modules\
`sort` Sorting the length of nucleotide sequences，the default parameters are shortest to longest,but --reverse option will reverses the result.\
`shuffle` Break up the order of the sequence.\
`blast` filter blast identify results. \
`extract`Randomly extracted sequences, seed number defaults to 10, with the option to enter your own seed (but it seems like I haven't gotten around to writing it yet)\
`n50` This option can read the fasta input and create a fai and output N50 and total length and GC content in a sample.txt\
`longest` can get the longest transcripts in cds sequence by gff3\
`RNAseq` can perform RNAseq analysis\


# **seqtools**
整合了一些用于阅读blast和hisat2的结果文件用于筛选和结果汇总的功能，并支持按geneid提取序列，为序列排序，随机提取或打乱序列顺序用于比对分析等功能。\
需求环境：python>3.5\
该工具依赖于Bio模块和pysam模块需要下载。\
对于linux 软件的依赖： Blast、samtools。 \
对于R包的依赖：limma、 edgeR、 DEseq2... \
查看帮助文档：`python main.py -h `\
目前seqtools有七个功能模块\
`sort` 对核苷酸序列的长度进行排序，默认参数为最短到最长，但--reverse选项将反转结果。\
`shuffle` 打乱序列顺序。\
`blast` 筛选blast比对率结果。\
`extract` 随机提取序列，种子数默认为10，但可以通过选项输入自己的种子（但我还没来得及写）。\
`N50` 此选项可以读取fasta输入并创建fai文件，输出N50、总长度和GC含量。\
`longest` 可以根据gff3文件获取cds序列的最长转录本。\
`RNAseq` 可以进行转录组分析的标准流程。 \

