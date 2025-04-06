import logging

import pysam

class MakeWindows:
    """
    建造窗口文件
    """
    def __init__(self, paths):
        self.paths = paths
    def build_fai(self):
        """
        Create FASTA index file (FAI).
        Output fai file for genome
        """
        try:
            logging.info("Creating FAI file...")
            with pysam.FastaFile(self.paths['fa_path']) as fasta:
                for chrom in fasta.references:
                    seq = fasta.fetch(chrom)
                    lines = seq.splitlines()
                    if len(lines) < 2:
                        continue
                    line_lengths = [len(line) for line in lines[:-1]]
                    if len(set(line_lengths)) > 1:
                        raise ValueError(f"Chromosome {chrom} has inconsistent line lengths: {line_lengths}")
            pysam.faidx(self.paths['fa_path'])
            logging.info("FAI file created successfully.")
        except Exception as e:
            logging.error(f"Error creating FAI file: {str(e)}")
            raise

    def calculate_chromosome_length(self):
        """
        Calculate and write chromosome lengths.This func read the fai and
            extract the first and second colunm to be a chromosome length file
        Input: .fai
        Output: chromosome length file
        """
        try:
            logging.info("Calculating chromosome lengths...")
            with open(self.paths['fai_path'], "r") as f:
                lines = [(line.strip().split("\t")[0], line.strip().split("\t")[1]) for line in f]
            with open(self.paths['chromosome_length'], "w") as f:
                for line_1, line_2 in lines:
                    f.write(f"{line_1}\t{line_2}\n")
            logging.info("Chromosome lengths written successfully.")
        except Exception as e:
            logging.error(f"Error calculating chromosome lengths: {str(e)}")
            raise

    def create_windows_py(self, windows_size=1000, step_size=500):
        """
        Create BED table of windows using chromosome lengths.
        input the windows size and step size for create
        """
        try:
            logging.info('Creating windows...')
            with open(self.paths['chromosome_length'], "r") as f:
                line = {lines.strip().split()[0]: lines.strip().split()[1] for lines in f}
                windows = []
                for key, length in line.items():
                    start = 0
                    while start < int(length):
                        windows.append((key, start, start + windows_size))
                        start += step_size
                with open(self.paths['windows_tsv'], 'w') as f:
                    for window in windows:
                        f.write(f'{window[0]}\t{window[1]}\t{window[2]}\n')
            logging.info('Windows created successfully.')
        except Exception as e:
            logging.error(f'Error creating windows: {str(e)}')
            raise

    def executor(self):
        # 执行函数
        self.build_fai()

        self.calculate_chromosome_length()

        self.create_windows_py()