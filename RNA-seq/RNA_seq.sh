#! /bin/bash

set -e  #error

if [ "$#" -ne 8 ]; then  #检测参数数量
   echo "
How to use this script
sh RNAseq_.sh <REF_genome> <gtf_file> <sample_file> <cpu> <method> <data_path> <script_path> <work_path>
 
   <REF_genome>    ref genome
   <gtf_file>      gtf file name
   <sample_file>   sample information table(sample.txt): the first colunm must is the prefix of sample which need to run RNA-seq
                   such as:
                   1st          2st
                   sample1       sample1.clean.fq.gz
                   sample2       sample2.clean.fq.gz
   <cpu>           hisat2 threads
   <method>        Variance analysis software， use edgeR or DESeq2
   <data_path>     RNA-seq data path（NOT ADD '/' IN THE TAIL）
   <script_path>   input your script path: like ~/path/to/script, the path need to have run-featurecounts.R abundance_estimates_to_matrix.pl and run_TMM_scale_matrix.pl
   <work_path>     result path
 Instruction manual：
 ##This is an important note, so please read it carefully##

 The RNA-seq data contain sample2.txt，the file need to contain:
 every path must use absolute path or relative path start by '~/' like '~/path/to/your/path'(but don't add '/' in the tail)
 Two columns, the first column is group class (custom) and the second column is the transcriptome data prefix:
   1st          2st
   group1        sample1.clean.fq.gz
   group1        sample2.clean.fq.gz
   group2        sample3.clean.fq.gz
   group2        sample4.clean.fq.gz
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
script=$7
work=$8
#检验参数是否正确
for THE_FILE in "$REF_GENOME" "$GTF_FILE" "$SAMPLE_FILE"; do
    if [ ! -f "$THE_FILE" ]; then
       echo "发生了一个错误，'$THE_FILE'文件不存在，或路径不正确。"
       exit 1

    fi
  done

#检验method参数是否正确
if [[ "$METHOD" != "edgeR" && "$METHOD" != "DESeq2" ]]; then
       echo "method参数设置错误，method参数请输入 edgeR 或 DESeq2 "
       exit 1
fi


#创建工作路径
mkdir -p $work/1.mapping $work/ref $work/2.Quantification $work/3.Merge_result $work/4.DE_analysis
#1.建立索引
cd $work/ref
if [ ! -f hisat2-build.log ]; then
	hisat2-build $REF_GENOME A  1>hisat2-build.log 2>&1 
fi

#2.转录组数据mapping到参考基因组上
cd $work/1.mapping
sample_count=$(awk 'NF{print $1}' $SAMPLE_FILE | uniq | wc -l)
sam_count=$(ls *.sam 2>/dev/null | wc -l)
if [ $sample_count -ne $sam_count ]; then
	awk '{print $1}' $SAMPLE_FILE | xargs -i echo "hisat2 --new-summary -p $CPU -x $work/ref/A -1 $DATA/{}_1.clean.fq.gz -2 $DATA/{}_2.clean.fq.gz -S {}.sam 1>{}.log 2>&1  " > run_hisat2.sh
       sh run_hisat2.sh	
fi
#3.sam文件转bam文件
cd $work/1.mapping
bam_count=$(ls *.bam 2>/dev/null | wc -l)
if [ $sample_count -ne $bam_count ]; then
awk '{print $1}' $SAMPLE_FILE | xargs -I {} sh -c 'samtools sort -o {}.bam {}.sam'

fi

#4.建立bam文件索引文件
cd $work/1.mapping
bai_count=$(ls *.bai 2>/dev/null | wc -l)
if [ $sample_count -ne $bai_count ]; then
awk '{print $1}' $SAMPLE_FILE | xargs -I {} sh -c 'samtools index {}.bam'

fi

#5.表达定量
cd $work/2.Quantification
if [ ! -f *.count ]; then
         awk '{print $1}' $SAMPLE_FILE | xargs -I {} sh -c "$script/run-featurecounts.R -b $work/1.mapping/{}.bam -g $GTF_FILE -o {} " > run_feature.sh
	 sh run_feature.sh
fi

#6.合并表达矩阵
cd $work/3.Merge_result
if [ ! -f *genes.TMM.EXPR.matrix ]; then
        ls $work/2.Quantification/*.count > genes.quant_files.txt
	perl $script/abundance_estimates_to_matrix.pl --est_method featureCounts --quant_files genes.quant_files.txt --out_prefix genes
fi

# 7.差异分析
cd $work/4.DE_analysis
perl /ME4012_Vol0002/bio_soft/public/miniconda3/bin/run_DE_analysis.pl --matrix $work/3.Merge_result/genes.counts.matrix --method $METHOD --samples_file $DATA/sample2.txt
