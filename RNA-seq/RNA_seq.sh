#! /bin/bash

set -e  #error

if [ "$#" -ne 6 ]; then  #检测参数数量
   echo "
How to use this script
sh RNAseq_.sh <REF_genome> <gtf_file> <sample_file> <cpu> <method> <data_path> 
 
   <REF_genome>    ref genome
   <gtf_file>      gtf file name
   <sample_file>   sample information table(sample.txt)
   <cpu>           hisat2 threads
   <method>        Variance analysis software， use edgeR or DESeq2
   <data_path>     RNA-seq data path（NOT ADD "/" IN THE TAIL）
 Instruction manual：
 ##This is an important note, so please read it carefully##
 The sample information table has the following format requirements：the first colunm contain RNA-seq sample's name。
 The RNA-seq data contain sample2.txt，the file need to contain:
 Two columns, the first column is group class (custom) and the second column is the transcriptome data prefix
  "

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
