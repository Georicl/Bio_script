#! /bin/bash
set -e #添加错误处理
if [ $# -ne 4 ] ;then
	echo "
脚本使用说明
sh sam_to_bam_delete.sh <sample> <index> <cpu> <RNA-seq>
sample	样本信息表，第一列为样本ID
index	ht2库的前缀路径
CPU	线程数
RNA_PATH RNA-seq	转录组数据路径（最后不要带/）
The purpose of this script is to automatically delete the SAM file after the alignment of a transcriptome file is completed, saving space.
	"
	exit 1
fi
# 导入参数

SAMPLE_FILE=$1
IDX=$2
CPU=$3
RNA_PATH=$4

# 遍历sample.txt提取的id
while read -r sample_name ;do
    #根据sample.txt定义提取的id
    FQ_1="${sample_name}_1.clean.fq.gz"
    FQ_2="${sample_name}_2.clean.fq.gz"
    SAM_FILE="${sample_name}.sam"
    BAM_FILE="${sample_name}.bam"
    #调用hisat2进行比对
    hisat2 --new-summary -x $IDX -p $CPU -1 "$4/$FQ_1" -2 "$4/$FQ_2" -S $SAM_FILE 1>"${sample_name}".log 2>&1
    #将SAM转换为BAM
    samtools sort -o $SAM_FILE $BAM_FILE
    rm $SAM_FILE

done < "$SAMPLE_FILE"
