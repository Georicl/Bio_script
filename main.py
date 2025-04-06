#!/usr/bin/env python3
"""
BIO-SeqTools: 生物信息学序列处理工具包
支持FASTA/BED文件操作、日志分析和质量评估
"""

import argparse
import sys
from pathlib import Path
from CQ_tools import man
from CQ_mapping import align
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

    # ==================== CQ工具命令 ==================== #
    @registry.register("cqmapping", "重测序数据比对")
    def _(parser):
        parser.formatter_class = argparse.RawTextHelpFormatter

        parser.description = """
        本程序使用了bwa进行重测序数据比对，使用samtools进行数据格式转换，如果用于进行后续CQ的分析。
        作者：向旸
        
        主要步骤包含：
        
        1. 比对 <- bwa
        2. 格式转换 <- samtools
        3. 构建索引 <- samtools 
        如果用于后续的CQ分析，根据雌雄测序的文件需要进行两次比对，我建议创建两个文件夹名为Female和Male，
        并在该目录下运行，以得到结果后直接使用cqtools进行比对。
        
        请不要在同一个输出目录的参数下运行两次该程序，该程序生成的结果文件名为：
        output.bam
        output.sort.bam
        output.sort.bam.bai
        因此请设立不同的输出目录进行比对计算。
        如果发生了意外的BUG，请在issue上提问，或直接发邮件至此处：Georicl@outlook.com
        If meet a bug when you run cqmapping, please ask in issue or email me:Georicl@outlook.com
        """
        parser.add_argument("--fasta", require=True, help="参考基因组输入路径")
        parser.add_argument("--pair1", require=True, help="双端测序文件1输入路径")
        parser.add_argument("--pair2", require=True, help="双端测序文件2输入路径")
        parser.add_argument("-t", "--threads", default=4, type=int, help="使用bwa进行比对的线程数量，默认为4")
        parser.add_argument("-o", "--output", required=True, help="输出路径目录")
        return lambda args: align.run_bwa_samtools(
            fasta=args.fasta,
            sample_1=args.pair1, sample_2=args.pair2,
            output_dir=args.output, cpu=args.threads
        )

    @registry.register("cqtools", "CQ值计算工具")
    def _(parser):
        parser.formatter_class = argparse.RawTextHelpFormatter

        parser.description = """
                本程序用于进行CQ的比对。
                作者：向旸

                主要步骤包含：

                1. 比对 <- bedtools
                2. 计算CQ
                3. 筛选
                如果使用cqmapping进行计算的情况下，可以使用如下参数进行计算：
                --f_bam Female/output.sort.bam
                --m_bam Male/output.sort.bam
                其他参数自行设置。
                最终结果文件为：F_M_CQ.filter.tsv
                
                如果发生了意外的BUG，请在issue上提问，或直接发邮件至此处：Georicl@outlook.com
                If meet a bug when you run cqmapping, please ask in issue or email me:Georicl@outlook.com
                """
        parser.add_argument("--f_bam", required=True, help="输入Female比对文件")
        parser.add_argument("--m_bam", required=True, help="输入Male比对文件")
        parser.add_argument("--fasta", required=True, help="输入参考基因组文件")
        parser.add_argument("--output", required=True, help="输出路径目录")
        parser.add_argument("--cq_value", type=float, default=0.3, help="CQ比对筛选阈值，默认为0.3，筛选范围为0-1")
        parser.add_argument("--parallel", type=int,choice=[1, 2], default=2, help="是否启用多进程进行bedtools覆盖度计算，选择值为1、2。默认值为2（启用），1为不启用")
        parser.add_argument("--reads_threshold", type=int, default=30, help="reads的支持数量，默认为30")
        return lambda args: man.main(
            f_bam=args.f_bam, m_bam=args.m_bam, fasta_path_get=args.fasta,
            output_path_get_path=args.output, cq_value=args.cq_value,
            num_parallel=args.parallel, reads_threshold=args.reads_threshold
        )



    # @registry.register("plot", "简易绘图")
    # def _(parser):
    #     parser.formatter_class = argparse.RawTextHelpFormatter
    #     parser.description = """\
    #             简易画图流程： 示例 main.py plot -R go/heatmap 参数
    #             如果是使用heatmap需要准备：
    #             genes.counts.matrix
    #             genes.TPM.not_cross_norm.TMM_info.txt
    #             *.DE_results
    #             genes.TMM.EXPR.matrix
    #             如果是使用go需要准备go富集或KEGG富集的通路excel表格：
    #             *.xlxs
    #             输出将直接输出为pdf，输出在output参数下选择的路径
    #             """
    #
    #     parser.add_argument("-R", "--Rs",type= str, help="选择进行何类画图,请选择输入go或heatmap，go画图可以输出go或kegg富集，heatmap输出热图和火山图")
    #     parser.add_argument("-s", "--sample", help="sample样本，如果是go请输入富集分析后得到的excel表格，如果是heatmap请输出genes.TPM.not_cross_norm.TMM_info.txt")
    #     parser.add_argument("-c", "--count", help="Heatmap选项：输入计数文件genes.counts.matrix")
    #     parser.add_argument("-e", "--expression", help="Heatmap选项：输入表达矩阵genes.TMM.EXPR.matrix")
    #     parser.add_argument("-d", "--deresult", help="Heatmap选项：*.DE_results")
    #     parser.add_argument("-m", "--map", default=None, help="Heatmap选项：对照的样本")
    #     parser.add_argument("-o", "--output", help="输出文件路径")
    #     return lambda args: seqtools.run_plot(
    #        args.R, args.s, args.count, args.expression, args.deresult, args.map, args.output
    #     )


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
