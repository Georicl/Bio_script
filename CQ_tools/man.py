import logging
import os

from CQ_tools import coverage_calculate, CQ_result, get_paths, merge_result, TMM, windows_make



def main(f_bam, m_bam, fasta_path_get, output_path_get_path, cq_value=0.3, num_parallel=2, reads_threshold=30):
    paths = get_paths.get_paths(f_bam=f_bam, m_bam=m_bam, fasta_path_get=fasta_path_get,
                                output_path_get_path=output_path_get_path)

    os.makedirs(paths['temp_dir'], exist_ok=True)

    logging.info("Step 1: Creating windows...")
    window_maker = windows_make.MakeWindows(paths)
    window_maker.executor()

    logging.info("Step 2: Splitting windows file...")
    temp_file_maker = coverage_calculate.MakeTempFile(paths)
    temp_file_maker.split_windws_file()

    logging.info("Step 3: Calculating coverage...")
    coverage_calculator = coverage_calculate.CoverageCalculate(paths, num_parallel=num_parallel)
    coverage_calculator.execute()

    logging.info("Step 4: Performing CQ calculation and filtering...")
    cq_calculator = CQ_result.CQCalCulate(paths, cq_value=cq_value, m_reads_threshold=reads_threshold)
    cq_calculator.cq_calculate()
    cq_calculator.cq_filter()
    cq_calculator.output_result()

    logging.info("All steps completed successfully!")

if __name__ == "__main__":
    import argparse
    import logging

    # 设置日志
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=[logging.StreamHandler()]
    )

    # 解析命令行参数
    parser = argparse.ArgumentParser(description="Run the complete CQ pipeline.")
    parser.add_argument("--f_bam", required=True, help="Path to the female BAM file.")
    parser.add_argument("--m_bam", required=True, help="Path to the male BAM file.")
    parser.add_argument("--fasta", required=True, help="Path to the reference FASTA file.")
    parser.add_argument("--output", required=True, help="Output directory path.")
    parser.add_argument("--cq_value", type=float, default=0.3, help="CQ filtering threshold (default: 0.3).")
    parser.add_argument("--parallel", type=int, default=4, help="Number of parallel processes (default: 4).")
    parser.add_argument("--reads_threshold", type=int, default=30, help="Number of reads threshold (default: 30).")
    args = parser.parse_args()
