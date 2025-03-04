#!/usr/bin/env python3
"""
BIO-SeqTools: 生物信息学序列处理工具包
支持FASTA/BED文件操作、日志分析和质量评估
"""

import argparse
import sys
from pathlib import Path

import seqtools


# region ################### 命令注册框架 #######################

class CommandRegistry:
    """命令行注册管理中心"""

    def __init__(self):
        self.handlers = {}
        self.parsers = {}
        self.root_parser = argparse.ArgumentParser(
            prog="bio_tools",
            description="生物信息学多功能处理工具",
            epilog="作者：向旸"
        )
        self.subparsers = self.root_parser.add_subparsers(
            title="可用命令",
            dest="command",
            required=True
        )

    def register(self, name, help_str):
        """命令注册装饰器工厂方法"""

        def decorator(config_func):
            parser = self.subparsers.add_parser(name, help=help_str)
            handler = config_func(parser)  # 配置参数并获取处理函数
            self.handlers[name] = handler
            self.parsers[name] = parser
            return handler

        return decorator


# endregion

# region ################### 具体命令实现 #######################
def setup_commands():
    """初始化所有子命令"""
    registry = CommandRegistry()

    # ======================== 核心功能命令 ======================== #

    @registry.register("sort", "按序列长度排序")
    def _(parser):
        """对应 seqtools.seq_sort"""
        parser.add_argument("gene_fasta", help="输入FASTA文件路径")
        parser.add_argument("out_path", help="输出目录路径")
        parser.add_argument("-r", "--reverse", action="store_true",
                            help="逆序排序（长到短）")
        return lambda args: seqtools.seq_sort(
            args.gene_fasta, args.out_path, args.reverse
        )

    @registry.register("shuffle", "随机打乱序列")
    def _(parser):
        """对应 seqtools.seq_shuffle"""
        parser.add_argument("gene_fasta", help="输入FASTA文件路径")
        parser.add_argument("out_path", help="输出目录路径")
        return lambda args: seqtools.seq_shuffle(
            args.gene_fasta, args.out_path
        )

    # ======================== 数据分析命令 ======================== #

    @registry.register("n50", "计算N50指标")
    def _(parser):
        """对应 seqtools.N50"""
        parser.add_argument("gene_fasta", help="输入FASTA文件路径")
        parser.add_argument("out_path", help="输出目录路径")
        parser.add_argument("-t", "--trinity", action="store_true",
                            help="Trinity格式输入文件")
        return lambda args: seqtools.N50(
            args.gene_fasta, args.out_path, args.trinity
        )

    @registry.register("blast", "过滤BLAST结果")
    def _(parser):
        """对应 seqtools.blast_identify"""
        parser.add_argument("blast_in", help="BLAST结果文件路径(format6)")
        parser.add_argument("out_path", help="输出目录路径")
        parser.add_argument("indentify", type=float,
                            help="相似度阈值（百分比单位）")
        return lambda args: seqtools.blast_identify(
            args.blast_in, args.out_path, args.indentify
        )

    @registry.register("longest", "提取最长转录本")
    def _(parser):
        parser.add_argument('gene_fasta', type=str, help='Input your genome **cds** FASTA file')
        parser.add_argument('out_path', type=str, help='Output path directory')
        parser.add_argument('gff', type=str, help='input gff3 file')
        return lambda args: seqtools.get_longest(
            args.gene_fasta, args.out_path, args.gff
        )

    # ==================== 其他实用工具命令 ==================== #

    @registry.register("extract", "随机抽取序列")
    def _(parser):
        """对应 seqtools.seq_extract"""
        parser.add_argument("gene_fasta", help="输入FASTA文件路径")
        parser.add_argument("out_path", help="输出目录路径")
        parser.add_argument("number", type=float,
                            help="抽取序列数量")
        parser.add_argument("-s", "--seed", type=int,
                            help="随机数种子")
        return lambda args: seqtools.seq_extract(
            args.gene_fasta, args.number, args.out_path, args.seed or 10
        )

    @registry.register("RNAseq", "转录组数据上游处理")
    def _(parser):
        # 使用 RawTextHelpFormatter 保留格式
        parser.formatter_class = argparse.RawTextHelpFormatter

        # 主描述（带格式的步骤说明）
        parser.description = """\
        处理RNA-seq数据分析完整流程，主要步骤包含：
        1. 参考基因组索引构建
        2. 序列比对（HISAT2）
        3. 表达量定量（featureCounts）
        4. 差异表达分析（DESeq2/edgeR）"""

        # 补充说明（带格式的文件示例）
        parser.epilog = """\
        【文件准备要求】

        > 样本信息表 (sample_path 参数文件):
        1st列: 样本名        （示例） | 2st列: 文件名 [可选]
        ------------------------------|-------------------------
        sample1               | sample1_1.clean.fq.gz
        sample1               | sample1_2.clean.fq.gz

        > 必须准备的 sample2.txt 文件:
        第一列: 组名          | 第二列: 样本名
        ------------------------------|-------------------------
        control               | sample1
        treatment             | sample2"""

        parser.add_argument('gene_fasta', help="输入参考基因组文件路径(以~/开头或绝对路径)")
        parser.add_argument('out_path', help="输出目录,该命令会在该目录下生成5个工作文件夹分别为ref，01-05")
        parser.add_argument('gtf_path', help="输入参考基因组gtf文件路径(以~/开头或绝对路径)")
        parser.add_argument('sample_path', help="输入样本信息表(以~/开头或绝对路径)，第一列为样本名（例如XYSM1_1.clean.fq.gz,样本名则为XYSM1，即_1/2...前的所有字符）")
        parser.add_argument('cpu', help="线程数")
        parser.add_argument('method', help="使用的定量方法，无生物学重复输入edgeR，有生物学重复输入DEseq2")
        parser.add_argument('data_path', help="输入转录组文件路径（末尾不带/）")
        return lambda args: seqtools.RNA_seq(
            Path(__file__).parent.absolute() / "RNA_seq" / "RNA_seq.sh",
            args.gene_fasta, args.gtf_path, args.sample_path, args.cpu, args.method, args.data_path,
            Path(__file__).parent.absolute() / "RNA_seq" / "support_script",
            args.out_path,
        )


    return registry





# endregion

# region #################### 主程序入口 #########################
def main():
    registry = setup_commands()

    try:
        # 空参数时显示帮助
        if len(sys.argv) == 1 :
            registry.root_parser.print_help()
            return

        if len(sys.argv) == 2 and sys.argv [1] in registry.parsers:
            registry.parsers[sys.argv[1]].print_help()  # 打印子命令帮助信息
            return

        # 参数解析与执行
        args = registry.root_parser.parse_args()
        if handler := registry.handlers.get(args.command):
            handler(args)
        else:
            raise RuntimeError(f"未知命令: {args.command}")

    except Exception as e:
        print(f"\033[31m[错误] {str(e)}\033[0m", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
# endregion
