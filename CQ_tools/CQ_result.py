import logging

import pandas as pd
from .TMM import process_tmm
import os


class CQCalCulate:
    def __init__(self, paths, cq_value=0.3, m_reads_threshold=30):
        self.paths = paths
        self.cq_value = float(cq_value)
        self.m_reads_threshold = int(m_reads_threshold)

    def cq_calculate(self):
        with open(self.paths['f_m_merge'], 'r') as f_in, \
                open(self.paths['merged_result'], 'w') as f_cov, \
                open(self.paths['reads_file'], 'w') as f_reads:  # 新增reads文件写入

            f_reads.write("window\tF\tM\n")  # 写入表头
            next(f_in)
            for lines in f_in:
                line = lines.strip().split('\t')
                chrom, start, end, f_reads_val, f_coverage, m_reads_val, m_coverage = line
                window = f"{chrom}-{start}-{end}"

                # 写入reads矩阵文件
                f_reads.write(f"{window}\t{f_reads_val}\t{m_reads_val}\n")

                # 原有计算逻辑
                m_coverage = float(m_coverage)
                f_coverage = float(f_coverage)
                cq_re = f_coverage / m_coverage if m_coverage != 0 else 0

                f_cov.write(
                    f"{chrom}\t{start}\t{end}\t{f_reads_val}\t{f_coverage}\t{m_reads_val}\t{m_coverage}\t{cq_re}\n")

    def _apply_normalization(self):
        """执行TMM标准化并重命名结果文件"""
        process_tmm(self.paths['reads_file'])  # 调用TMM模块
        src = f"{self.paths['reads_file']}_TMM.txt"  # 生成的默认路径
        os.rename(src, self.paths['tmm_normalized'])  # 重命名到目标路径

    def cq_filter(self):
        # 执行标准化并校验数据
        self._apply_normalization()

        try:
            # 读取标准化数据（保留窗口列）
            normalized_df = pd.read_csv(self.paths['tmm_normalized'], sep='\t')

            # 读取原始合并结果并创建窗口ID
            merged_df = pd.read_csv(
                self.paths['merged_result'],
                sep='\t',
                names=['chrom', 'start', 'end', 'f_reads', 'f_cov', 'm_reads', 'm_cov', 'cq_re']
            )
            merged_df['window'] = merged_df.apply(lambda x: f"{x['chrom']}-{x['start']}-{x['end']}", axis=1)

            # 数据合并校验
            original_length = len(merged_df)
            combined_df = pd.merge(
                merged_df,
                normalized_df[['window', 'M']],
                on='window',
                how='inner'
            )

            # 记录丢失窗口数量（容错日志）
            if len(combined_df) < original_length:
                lost = original_length - len(combined_df)
                print(f"警告：标准化数据缺失 {lost} 个窗口，已自动过滤")

            # 应用过滤条件
            filtered_df = combined_df[
                (combined_df.cq_re <= self.cq_value) &
                (combined_df.M >= 30)
                ]

            # 替换标准化值并保存
            filtered_df['m_reads'] = filtered_df['M']  # 直接替换原m_reads列为标准化值
            filtered_df.to_csv(
                self.paths['cq_output_f'],
                sep='\t',
                columns=['chrom', 'start', 'end', 'f_reads', 'f_cov', 'm_reads', 'm_cov', 'cq_re'],
                index=False,
                header=False
            )

        except Exception as e:
            print(f"处理失败：{str(e)}")
            raise

    def output_result(self):
        """筛选cq_output_f文件，输出最终结果到filtered_cq"""
        try:
            # 读取cq_output_f文件
            df = pd.read_csv(
                self.paths['cq_output_f'],
                sep='\t',
                names=['chrom', 'start', 'end', 'f_reads', 'f_cov', 'm_reads', 'm_cov', 'cq_re']
            )

            # 应用过滤条件
            filtered_df = df[
                (df['cq_re'] < self.cq_value) &
                (df['m_reads'] > self.m_reads_threshold)
                ]

            # 保存最终结果
            filtered_df.to_csv(
                self.paths['filtered_cq'],
                sep='\t',
                columns=['chrom', 'start', 'end', 'f_reads', 'f_cov', 'm_reads', 'm_cov', 'cq_re'],
                index=False,
                header=False
            )

            logging.info(f"最终结果已保存到 {self.paths['filtered_cq']}")

        except Exception as e:
            logging.error(f"输出结果失败: {str(e)}")

