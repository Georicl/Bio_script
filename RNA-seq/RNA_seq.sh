#! /bin/bash

set -e  #error

if [ "$#" -ne 6 ]; then  #检测参数数量
   echo "
脚本使用方法
sh RNAseq_.sh <REF_genome> <gtf_file> <sample_file> <cpu> <method> <data_path> 
 
   <REF_genome>    参考基因组文件名
   <gtf_file>      注释文件（gtf）名
   <sample_file>   样本信息表（如sample.txt）文件名
   <cpu>           hisat2使用的线程数
   <method>        差异分析所使用的软件，请输入edgeR（无生物学重复）或DESeq2（有生物学重复）
   <data_path>     转录组数据存放路径
 使用说明：
 ##该说明十分重要，请将该说明仔细看完## 
 样本信息表进行如下格式要求： 第一列必须为转录组数据_1/2.clean.fq.gz的前缀，即转录组数据名。
 1.在转录组数据中需要包含一个sample2.txt，该文件应包含如下内容：
为两列，第一列为组类（自定义），第二列为转录组数据前缀名
 2.转录组数据应至少包括以下文件：
信息表文件 sample2.txt 参考基因组序列文件 参考基因注释文件 原始测序数据文件"
   exit 1
fi

#变量声明
REF_GENOME=$1
GTF_FILE=$2
SAMPLE_FILE=$3
CPU=$4
METHOD=$5
DATA=$6    
#检验参数是否正确
for THE_FILE in "$REF_GENOME" "$GTF_FILE" "$SAMPLE_FILE"; do
    if [ ! -f "$THE_FILE" ]; then
       echo "发生了一个错误，'$THE_FILE'文件不存在，或路径不正确，请将文件放置在该脚本同一文件夹下。"
       exit 1

    fi
  done

#检验method参数是否正确
if [[ "$METHOD" != "edgeR" && "$METHOD" != "DESeq2" ]]; then
       echo "method参数设置错误，method参数请输入 edgeR 或 DESeq2 "
       exit 1
fi


#创建工作路径
mkdir -p ~/Yang/work/1.mapping ~/Yang/work/ref ~/Yang/work/2.Quantification ~/Yang/work/3.Merge_result ~/Yang/work/4.DE_analysis
#1.建立索引
cd ~/Yang/work/ref
if [ ! -f "hisat2-build.log" ]; then
	hisat2-build $REF_GENOME A  1>hisat2-build.log 2>&1 
fi

#2.转录组数据mapping到参考基因组上
cd ~/Yang/work/1.mapping
if [ ! -f "*.sam" ]; then
	awk '{print $2}' $SAMPLE_FILE | xargs -i echo "hisat2 --new-summary -p $CPU -x ~/Yang/work/ref/A -1 $DATA/{}_1.clean.fq.gz -2 $DATA/{}_2.clean.fq.gz -S {}.sam 1>{}.log 2>&1  " > run_hisat2.sh
       sh run_hisat2.sh	
fi
#3.sam文件转bam文件
cd ~/Yang/work/1.mapping
if [ ! -f "*.bam" ]; then      
awk '{print $2}' $SAMPLE_FILE | xargs -I {} sh -c 'samtools sort -o {}.bam {}.sam' 

fi

#4.建立bam文件索引文件
cd ~/Yang/work/1.mapping
if [ ! -f "*.bai" ]; then
awk '{print $2}' $SAMPLE_FILE | xargs -I {} sh -c 'samtools index {}.bam'

fi

#5.表达定量
cd ~/Yang/work/2.Quantification
if [ ! -f "*.count" ]; then
         awk '{print $2}' $SAMPLE_FILE | xargs -I {} sh -c "~/script/run-featurecounts.R -b ~/Yang/work/1.mapping/{}.bam -g $GTF_FILE -o {} " > run_feature.sh
	 sh run_feature.sh
fi

#6.合并表达矩阵
cd ~/Yang/work/3.Merge_result
if [ ! -f "*genes.TMM.EXPR.matrix" ]; then
        ls ~/Yang/work/2.Quantification/*.count > genes.quant_files.txt
	perl ~/script/abundance_estimates_to_matrix.pl --est_method featureCounts --quant_files genes.quant_files.txt --out_prefix genes
fi

# 7.差异分析
cd ~/Yang/work/4.DE_analysis
perl /ME4012_Vol0002/bio_soft/public/miniconda3/bin/run_DE_analysis.pl --matrix ~/Yang/work/3.Merge_result/genes.counts.matrix --method $METHOD --samples_file $DATA/sample2.txt
