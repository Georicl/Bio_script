import pysam
import matplotlib.pyplot as plt

# ========== 设置输入文件 ==========
# BAM 文件及索引文件
bam_file = "path/to/your_file.bam"  # 替换为你的 BAM 文件路径
index_file = "path/to/your_file.bam.bai"  # 替换为 BAM 的索引文件路径

# 设置染色体、起始和终止位置
chromosome = "chr1"
start_pos = 100000
end_pos = 200000


# ========== 加载 BAM 文件并绘制覆盖度 ==========
def plot_bam_coverage(bam_file, chromosome, start_pos, end_pos):
    """使用 pysam 绘制 BAM 文件覆盖度"""
    samfile = pysam.AlignmentFile(bam_file, "rb")
    coverage = []
    positions = []

    # 遍历指定区域的覆盖度
    for pileup_column in samfile.pileup(chromosome, start_pos, end_pos):
        positions.append(pileup_column.reference_pos)
        coverage.append(pileup_column.n)

    # 绘制覆盖度曲线
    plt.figure(figsize=(10, 6))
    plt.plot(positions, coverage, label="BAM Coverage", color="blue")
    plt.xlabel("Genomic Position")
    plt.ylabel("Coverage")
    plt.title(f"BAM Coverage for {chromosome}:{start_pos}-{end_pos}")
    plt.legend()
    plt.grid()
    plt.show()


# ========== 主函数 ==========
if __name__ == "__main__":
    plot_bam_coverage(bam_file, chromosome, start_pos, end_pos)