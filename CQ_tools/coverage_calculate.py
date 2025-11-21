import logging
import os
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor
from concurrent.futures.process import ProcessPoolExecutor


class MakeTempFile:
    # 分割临时文件
    def __init__(self, paths):
        self.paths = paths

    def split_windws_file(self):
        with open(self.paths['windows_tsv'], 'r') as f_in:
            windows = f_in.readlines()

        chrom_windows = {}
        for window in windows:
            chrom, start, end = window.strip().split('\t')
            if chrom not in chrom_windows:
                chrom_windows[chrom] = []
            chrom_windows[chrom].append((start, end))

        temp_dir = self.paths['temp_dir']
        if not os.path.exists(temp_dir):
            os.makedirs(temp_dir)

        for chrom, windows in chrom_windows.items():
            file_path = os.path.join(temp_dir, f"{chrom}.bed")
            with open(file_path, "w") as f_out:
                for start, end in windows:
                    f_out.write(f"{chrom}\t{start}\t{end}\n")
            logging.info(f"Create file for {chrom}: {file_path}")

class CoverageCalculate:
    def __init__(self, paths, num_parallel=4):
        self.paths = paths
        self.parallel = int(num_parallel)

    def _run_full_coverage(self, bam_input, output_path):
        """运行 bedtools coverage 对原始窗口文件计算覆盖度，并直接写入最终路径"""
        cmd = (
            f"bedtools coverage -a {self.paths['windows_tsv']} -b {bam_input} -sorted "
            f"> {output_path}"
        )
        try:
            subprocess.run(cmd, shell=True, check=True)
            logging.info(f"完整覆盖度文件生成完成: {output_path}")
        except subprocess.CalledProcessError as e:
            logging.error(f"完整覆盖度文件生成失败: {e}")

    def m_coverage_calculate(self) -> None:
        """优化后的 M-bam 覆盖度计算"""
        start_time = time.time()

        # 直接生成到最终路径
        self._run_full_coverage(self.paths['m_bam'], self.paths['m_cov'])

        logging.info(f"M-bam 计算完成，耗时: {time.time() - start_time:.2f}秒")

    def f_coverage_calculate(self) -> None:
        """优化后的 F-bam 覆盖度计算"""
        start_time = time.time()

        # 直接生成到最终路径
        self._run_full_coverage(self.paths['f_bam'], self.paths['f_cov'])

        logging.info(f"F-bam 计算完成，耗时: {time.time() - start_time:.2f}秒")

    def _merge_final_result(self) -> None:
        """合并雌雄结果并校验每个窗口的chr/start/end完全一致"""
        header = "chr\tstart\tend\tf_num_reads\tf_coverage\tm_num_reads\tm_coverage\n"
        mismatch_count = 0
        line_count = 0

        with open(self.paths['f_m_merge'], "w") as f_out:
            f_out.write(header)

            # 一次性读取文件内容
            with open(self.paths['f_cov'], "r") as f_f, \
                    open(self.paths['m_cov'], "r") as f_m:
                lines_f = f_f.readlines()
                lines_m = f_m.readlines()

            # 跳过表头
            lines_f = lines_f[1:]
            lines_m = lines_m[1:]

            # 遍历并合并内容
            for line_num, (line_f, line_m) in enumerate(zip(lines_f, lines_m), start=1):
                parts_f = line_f.strip().split("\t")
                parts_m = line_m.strip().split("\t")

                # 提取区域信息
                chr_f, start_f, end_f = parts_f[0], int(parts_f[1]), int(parts_f[2])
                chr_m, start_m, end_m = parts_m[0], int(parts_m[1]), int(parts_m[2])

                # 严格校验区域一致性
                if chr_f != chr_m or start_f != start_m or end_f != end_m:
                    logging.error(
                        f"行 {line_num} 区域不匹配:\n"
                        f"雌性文件: {chr_f}:{start_f}-{end_f}\n"
                        f"雄性文件: {chr_m}:{start_m}-{end_m}"
                    )
                    mismatch_count += 1
                    continue  # 跳过不匹配行

                # 提取覆盖度数据（假设bedtools输出列顺序固定）
                f_num_reads = parts_f[3]
                f_coverage = parts_f[6]
                m_num_reads = parts_m[3]
                m_coverage = parts_m[6]

                # 写入合并行
                merged_line = (
                    f"{chr_f}\t{start_f}\t{end_f}\t"
                    f"{f_num_reads}\t{f_coverage}\t"
                    f"{m_num_reads}\t{m_coverage}\n"
                )
                f_out.write(merged_line)
                line_count += 1

            # 统计剩余行数
            remaining_f = len(lines_f[line_count:])
            remaining_m = len(lines_m[line_count:])

            if remaining_f > 0 or remaining_m > 0:
                logging.error(
                    f"文件行数不一致！雌性文件剩余 {remaining_f} 行，雄性文件剩余 {remaining_m} 行"
                )
            if mismatch_count > 0 or remaining_f + remaining_m > 0:
                logging.warning(f"合并完成，但发现 {mismatch_count} 处区域不匹配和 {remaining_f + remaining_m} 行缺失")
            else:
                logging.info(f"成功合并 {line_count} 行数据至 {self.paths['f_m_merge']}")

    def execute(self) -> None:
        """
        计算覆盖度，雌雄一起计算
        """
        start_time = time.time()

        with ProcessPoolExecutor(max_workers=min(2, self.parallel)) as executor:
            future_m = executor.submit(self.m_coverage_calculate)
            future_f = executor.submit(self.f_coverage_calculate)

            # 等待两个任务完成
            future_m.result()
            future_f.result()

        self._merge_final_result()

        logging.info(f"覆盖度计算完成，总耗时: {time.time() - start_time:.2f}秒")

# class CoverageCalculate:
#     # 计算主函数
#     def __init__(self, paths, num_parallel=4):
#         self.paths = paths
#         self.parallel = int(num_parallel)
#     def _run_full_coverage(self, bam_input, output_path):
#         """运行 bedtools coverage 对原始窗口文件计算覆盖度"""
#         cmd = (
#             f"bedtools coverage -a {self.paths['windows_tsv']} -b {bam_input} -sorted "
#             f"> {output_path}"
#         )
#         try:
#             subprocess.run(cmd, shell=True, check=True)
#             logging.info(f"完整覆盖度文件生成完成: {output_path}")
#         except subprocess.CalledProcessError as e:
#             logging.error(f"完整覆盖度文件生成失败: {e}")
#
#     def _split_coverage_by_chrom(self, full_coverage_path, output_suffix):
#         """按染色体分割完整覆盖度文件"""
#         chrom_files = {}
#         with open(full_coverage_path, "r") as f_in:
#             for line in f_in:
#                 parts = line.strip().split("\t")
#                 chrom = parts[0]
#                 if chrom not in chrom_files:
#                     chrom_files[chrom] = open(
#                         os.path.join(self.paths['temp_dir'], f"{chrom}_coverage_{output_suffix}.bed"), "w"
#                     )
#                 chrom_files[chrom].write(line)
#
#         # 关闭所有文件句柄
#         for chrom, f_out in chrom_files.items():
#             f_out.close()
#             logging.info(f"按染色体分割完成: {chrom}")
#
#
#     def _get_chrom_list_sorted(self):
#         """获取按窗口数升序排列的染色体列表（小染色体优先）"""
#         chrom_sizes = {}
#         for f in os.listdir(self.paths['temp_dir']):
#             if f.endswith(".bed") and not f.endswith(("_coverage_f.bed", "_coverage_m.bed")):
#                 chrom = os.path.splitext(f)[0]
#                 file_path = os.path.join(self.paths['temp_dir'], f)
#                 with open(file_path, "r") as bed_file:
#                     line_count = sum(1 for _ in bed_file)
#                 chrom_sizes[chrom] = line_count
#
#         # 按窗口数升序排序（小染色体在前）
#         sorted_chroms = sorted(chrom_sizes.items(), key=lambda x: x[1])
#         return [chrom for chrom, _ in sorted_chroms]
#
#     def m_coverage_calculate(self) -> None:
#         """优化后的 M-bam 覆盖度计算"""
#         start_time = time.time()
#
#         # 生成完整覆盖度文件
#         full_coverage_path = os.path.join(self.paths['temp_dir'], "full_coverage_m.bed")
#         self._run_full_coverage(self.paths['m_bam'], full_coverage_path)
#
#         # 按染色体分割覆盖度文件
#         self._split_coverage_by_chrom(full_coverage_path, output_suffix="m")
#
#         logging.info(f"M-bam 计算完成，耗时: {time.time() - start_time:.2f}秒")
#
#     def f_coverage_calculate(self) -> None:
#         """优化后的 F-bam 覆盖度计算"""
#         start_time = time.time()
#
#         # 生成完整覆盖度文件
#         full_coverage_path = os.path.join(self.paths['temp_dir'], "full_coverage_f.bed")
#         self._run_full_coverage(self.paths['f_bam'], full_coverage_path)
#
#         # 按染色体分割覆盖度文件
#         self._split_coverage_by_chrom(full_coverage_path, output_suffix="f")
#
#         logging.info(f"F-bam 计算完成，耗时: {time.time() - start_time:.2f}秒")
#
#     def execute(self) -> None:
#         """
#         计算覆盖度，雌雄一起计算
#         """
#         self.m_coverage_calculate()
#         self.f_coverage_calculate()