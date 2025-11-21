import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# ============================ 配置区域 ============================
# 输入文件路径 (请修改为您实际的文件名)
rm_file = '/Volumes/biomatics/GDOU-work/chr4_col/old_te_ann/SanyaOcu.p_ctg.FINAL.fasta.out'

# 目标染色体名称 (请确保与文件第一列一致，如 'Chr4', 'sy.chr4', '4' 等)
target_chrom_id = 'Chr4'

# 感兴趣的 Gap 区域 (SDR Candidate)
gap_start = 1
gap_end = 6800000

# 滑动窗口设置 (用于密度图)
window_size = 50000  # 50kb 窗口
step_size = 10000  # 10kb 步长

# 绘图范围 (在 Gap 前后各多画多少 bp)
plot_buffer = 2000000


# ================================================================

def parse_repeatmasker_out(filename):
    """
    解析 RepeatMasker 标准 .out 文件
    """
    print(f"[INFO] 正在读取文件: {filename} ...")
    try:
        # 读取 .out 文件，跳过前3行头信息
        # 列名根据标准 RepeatMasker 输出定义
        df = pd.read_csv(filename, sep=r'\s+', skiprows=3, header=None, engine='python',
                         names=['SW_score', 'perc_div', 'perc_del', 'perc_ins', 'query_seq',
                                'pos_start', 'pos_end', 'pos_left', 'strand',
                                'matching_repeat', 'repeat_class_family',
                                'r_start', 'r_end', 'r_left', 'ID', 'overlap_flag'])

        # 简单的数据清洗
        df['perc_div'] = pd.to_numeric(df['perc_div'], errors='coerce')
        df['pos_start'] = pd.to_numeric(df['pos_start'], errors='coerce')
        df['pos_end'] = pd.to_numeric(df['pos_end'], errors='coerce')

        return df
    except Exception as e:
        print(f"[ERROR] 读取失败: {e}")
        return None


def classify_te(family_str):
    """
    将细分的重复序列家族归类为四大类
    """
    fam = str(family_str).upper()
    if 'UNK' in fam or 'UNCLASSIFIED' in fam:  # 捕捉 Unclassified
        return 'Unclassified'
    if 'LTR' in fam or 'GYPSY' in fam or 'COPIA' in fam:
        return 'LTR'
    elif 'LINE' in fam:
        return 'LINE'
    elif 'DNA' in fam or 'HELITRON' in fam:
        return 'DNA'
    else:
        return 'Other'


# ============================ 主程序 ============================

# 1. 数据加载
df = parse_repeatmasker_out(rm_file)

if df is not None:
    # 2. 筛选目标染色体
    # 使用模糊匹配以防名称差异
    df_chrom = df[df['query_seq'].astype(str).str.contains(target_chrom_id, case=False)].copy()

    if df_chrom.empty:
        print(f"[WARNING] 未找到染色体 '{target_chrom_id}'。请检查文件中的染色体命名。")
        print("文件中的染色体名示例:", df['query_seq'].unique()[:5])
    else:
        print(f"[INFO] 筛选出 {len(df_chrom)} 条位于 {target_chrom_id} 的记录。")

        # 添加大类分类
        df_chrom['SuperFamily'] = df_chrom['repeat_class_family'].apply(classify_te)

        # ================= 图表 1: TE 密度景观图 =================
        print("[INFO] 正在生成 TE 密度数据...")

        plot_start = max(0, gap_start - plot_buffer)
        plot_end = gap_end + plot_buffer

        # 仅分析绘图范围内的数据以加速
        df_roi = df_chrom[(df_chrom['pos_end'] >= plot_start) & (df_chrom['pos_start'] <= plot_end)]

        landscape_data = []
        for pos in range(plot_start, plot_end, step_size):
            w_s = pos
            w_e = pos + window_size

            # 简单统计：找出中心点在窗口内的 TE
            centers = (df_roi['pos_start'] + df_roi['pos_end']) / 2
            in_window = df_roi[(centers >= w_s) & (centers < w_e)]

            # 计算各类型覆盖长度
            total_len = (in_window['pos_end'] - in_window['pos_start']).sum()
            ltr_len = (in_window[in_window['SuperFamily'] == 'LTR']['pos_end'] -
                       in_window[in_window['SuperFamily'] == 'LTR']['pos_start']).sum()
            line_len = (in_window[in_window['SuperFamily'] == 'LINE']['pos_end'] -
                        in_window[in_window['SuperFamily'] == 'LINE']['pos_start']).sum()
            unclass_len = (in_window[in_window['SuperFamily'] == 'Unclassified']['pos_end'] -
                           in_window[in_window['SuperFamily'] == 'Unclassified']['pos_start']).sum()

            landscape_data.append({
                'Position_Mb': w_s / 1e6,
                'Total_TE': min(1.0, total_len / window_size),
                'LTR': min(1.0, ltr_len / window_size),
                'LINE': min(1.0, line_len / window_size),
                'Unclassified': min(1.0, unclass_len / window_size)
            })

        df_plot = pd.DataFrame(landscape_data)

        plt.figure(figsize=(14, 6))
        plt.plot(df_plot['Position_Mb'], df_plot['Total_TE'], color='lightgrey', label='Total Repeats', alpha=0.5)
        plt.plot(df_plot['Position_Mb'], df_plot['LTR'], color='#e74c3c', label='LTR Retrotransposons', linewidth=2)
        plt.plot(df_plot['Position_Mb'], df_plot['LINE'], color='#3498db', label='LINEs', linewidth=1.5, alpha=0.7)
        plt.plot(df_plot['Position_Mb'], df_plot['Unclassified'], color='purple', label='Unclassified Repeats',
                 linewidth=2)
        # 标记 Gap
        plt.axvspan(gap_start / 1e6, gap_end / 1e6, color='#f1c40f', alpha=0.3, label='Synteny Gap (Target Region)')

        plt.title(f'TE Density Landscape: {target_chrom_id} ({plot_start // 1000}k - {plot_end // 1000}k)', fontsize=14)
        plt.xlabel('Position (Mb)', fontsize=12)
        plt.ylabel('Coverage Fraction', fontsize=12)
        plt.legend(loc='upper right')
        plt.grid(True, linestyle=':', alpha=0.6)
        plt.ylim(0, 1.0)
        plt.tight_layout()
        plt.show()

        # ================= 图表 2: LTR 插入时间 (差异度) 分析 =================
        print("[INFO] 正在对比 Gap 区域内外的 LTR 差异度 (插入时间)...")

        # 定义区域
        # Inside Gap: 位于 Gap 内部的 TE
        mask_inside = (df_chrom['pos_start'] >= gap_start) & (df_chrom['pos_end'] <= gap_end)

        # Outside Gap: 位于 Gap 外部（作为对照背景）
        # 这里取整条染色体的其他部分，或者您可以限制在 Gap 附近的 Flanking 区域
        mask_outside = ~mask_inside

        ltr_inside = df_chrom[mask_inside & (df_chrom['SuperFamily'] == 'LTR')]
        ltr_outside = df_chrom[mask_outside & (df_chrom['SuperFamily'] == 'LTR')]

        if len(ltr_inside) > 0:
            plt.figure(figsize=(10, 6))

            # 绘制直方图 (KDE + Hist)
            sns.histplot(ltr_outside['perc_div'], color='grey', label='Background (Outside Gap)',
                         kde=True, stat="density", element="step", alpha=0.3, bins=50)
            sns.histplot(ltr_inside['perc_div'], color='red', label='Target Region (Inside Gap)',
                         kde=True, stat="density", element="step", alpha=0.6, bins=50)

            plt.title('LTR Divergence Distribution (Insertion Age Proxy)', fontsize=14)
            plt.xlabel('Sequence Divergence (%) \n(Low Divergence = Recent Insertion / Young)', fontsize=12)
            plt.ylabel('Density', fontsize=12)
            plt.legend()
            plt.grid(True, alpha=0.3)
            plt.xlim(0, 40)  # 通常关注 0-40% 的差异度

            plt.tight_layout()
            plt.show()

            print("\n[结果解读提示]")
            print("图表 2 (差异度分布):")
            print(" - 如果红色曲线的峰值显著偏左 (靠近 0-5% 差异度)，说明该区域有近期爆发的年轻 LTR。")
            print("   -> 结论: 这是一个新形成的、正在快速演化的性染色体区域。")
            print(" - 如果红色曲线与灰色曲线重叠或偏右，说明该区域是古老的重复序列堆积。")
        else:
            print("[INFO] Gap 区域内没有检测到 LTR，无法进行差异度对比。")