import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import pandas as pd
import numpy as np
regions_df = pd.read_csv('chr4.txt', sep='\t', header=None, names=['chromosome', 'start', 'end', 'annotation'])

### 全局字体 (Global Font Settings)
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial']

### 数据定义 (Data Definitions)

## chromosome
chr4_len = 45703401

## amhy
amhy_start = 314775
amhy_end = 317071
# 扩增产物 (Amplicon)
amhy_primer_f_relative_start = 578
amhy_product_len = 662
amhy_primer_start = amhy_start + amhy_primer_f_relative_start - 1
amhy_primer_end = amhy_primer_start + amhy_product_len - 1

## amh
amh_start = 19851471
amh_end = 19853845
# 扩增产物 (Amplicon)
amh_primer_f_relative_start = 605
amh_product_len = 742
amh_primer_start = amh_start + amh_primer_f_relative_start - 1
amh_primer_end = amh_primer_start + amh_product_len - 1


### 图像绘制 (Plotting)
# 底图 (Base Figure)
fig, axes = plt.subplots(
    4, 1,
    figsize=(15, 10),
    gridspec_kw={'height_ratios': [1, 1, 1,4]}
)

### amhy
ax_amhy = axes[0]
ax_amhy.set_title('amhy', fontsize=14, loc='left', style='italic')
ax_amhy.hlines(y=0, xmin=amhy_start, xmax=amhy_end, color='#ad002b', linewidth=8, alpha=0.8, label='amhy')
ax_amhy.hlines(y=0, xmin=amhy_primer_start, xmax=amhy_primer_end, colors='gray', linewidth=8, label=f'Amplicon ({amhy_product_len} bp)')

# ==================== MODIFICATION START ====================
# 计算amhy的总长度 (Calculate total length of amhy)
amhy_len = amhy_end - amhy_start
# 在基因条的右侧添加总长度文本 (Add total length text to the right of the gene bar)
ax_amhy.text(
    amhy_end + 100, 0,  # Position: a little right of the end, vertically centered
    f'{amhy_len} bp',   # Text string
    ha='left',          # Horizontal alignment
    va='center',        # Vertical alignment
    fontsize=12
)
# ===================== MODIFICATION END =====================

arrow_props = dict(facecolor='black', shrink=0.05, width=1, headwidth=8)
ax_amhy.annotate(
    '', xy=(amhy_primer_start, 0.05),
    xytext=(amhy_primer_start, 0.3),
    arrowprops=arrow_props
)
ax_amhy.annotate(
    '', xy=(amhy_primer_end, -0.05),
    xytext=(amhy_primer_end, -0.3),
    arrowprops=arrow_props
)
ax_amhy.text(
    amhy_primer_start - 600, 0.4,
    'Sysex-F: GCCGCACTCACAGGTGAAGA',
    ha='left', va='center', fontsize=11, rotation=0
)
ax_amhy.text(
    amhy_primer_end + 600, - 0.4,
    'Sysex-R: CCACCAATGATGGCGTCTAT',
    ha='right', va='center', fontsize=11, rotation=0
)

ax_amhy.text(
    (amhy_primer_start + amhy_primer_end) / 2, 0.1,
    f'{amhy_product_len} bp',
    ha='center', va='bottom', fontsize=15, fontweight='bold'
)

padding = (amhy_end - amhy_start) * 0.1
ax_amhy.set_xlim(amhy_start - padding, amhy_end + padding)
ax_amhy.set_ylim(-1, 1)
ax_amhy.ticklabel_format(useOffset=False, style='plain')
ax_amhy.set_yticks([])
ax_amhy.spines['top'].set_visible(False)
ax_amhy.spines['right'].set_visible(False)
ax_amhy.spines['left'].set_visible(False)


### 染色体图 (Chromosome Map)
ax_chr = axes[1]

# 染色体背景(乳白色) (Chromosome background)
ax_chr.hlines(y=0, xmin=1, xmax= chr4_len, color='#f5f3f0', linewidth=10)
# CQ信息(灰色) (CQ information)
for index, row in regions_df.iterrows():
    ax_chr.hlines(y=0, xmin=row['start'], xmax=row['end'], color='gray', linewidth=10)

# 标记amhy和amh坐标 (Mark amhy and amh coordinates)
ax_chr.vlines(x=[amhy_start, amh_start], ymin=-0.2, ymax=0.2, colors=['#ad002b', '#00468c'], linewidth=3, clip_on=False)
# amh/amhy标签 (amh/amhy labels)
ax_chr.text(amhy_start, 0.3, 'amhy', ha='center', va='bottom', fontsize=12, color='#ad002b', style='italic')
ax_chr.text(amh_start, 0.3, 'amh', ha='center', va='bottom', fontsize=12, color='#00468c', style='italic')

# 坐标轴 (Axes)
ax_chr.set_xlim(0, chr4_len)
ax_chr.set_ylim(-0.03, 1.2)
# 修改单位坐标 (Modify coordinate units)
ax_chr.set_xlabel('Chromosome Position (kb)')
formatter = mticker.FuncFormatter(lambda x, p: f'{int(x/1000):,}')
ax_chr.xaxis.set_major_formatter(formatter)

ticks = ax_chr.get_xticks()
ticks = ticks[ticks <= chr4_len]
ticks = np.append(ticks, chr4_len)
ticks = np.unique(ticks)
ax_chr.set_xticks(ticks)
ax_chr.set_yticks([])

# 边框设置 (Border settings)
ax_chr.spines['top'].set_visible(False)
ax_chr.spines['right'].set_visible(False)
ax_chr.spines['left'].set_visible(False)

### amh
ax_amh = axes[2]
ax_amh.set_title('amh', fontsize=14, loc='left', style='italic')
ax_amh.hlines(y=0, xmin=amh_start, xmax=amh_end, color='#00468c', linewidth=8, alpha=0.8, label='amh')
ax_amh.hlines(y=0, xmin=amh_primer_start, xmax=amh_primer_end, colors='gray', linewidth=8, label=f'Amplicon ({amh_product_len} bp)')

# ==================== MODIFICATION START ====================
# 计算amh的总长度 (Calculate total length of amh)
amh_len = amh_end - amh_start
# 在基因条的右侧添加总长度文本 (Add total length text to the right of the gene bar)
ax_amh.text(
    amh_end + 100, 0,    # Position: a little right of the end, vertically centered
    f'{amh_len} bp',     # Text string
    ha='left',           # Horizontal alignment
    va='center',         # Vertical alignment
    fontsize=12
)
# ===================== MODIFICATION END =====================

ax_amh.annotate(
    '', xy=(amh_primer_start, 0.05),
    xytext=(amh_primer_start, 0.3),
    arrowprops=arrow_props
)
ax_amh.annotate(
    '', xy=(amh_primer_end, -0.05),
    xytext=(amh_primer_end, -0.3),
    arrowprops=arrow_props
)
ax_amh.text(
    amh_primer_start - 600, 0.4,
    'Sysex-F: GCCGCACTCACAGGTGAAGA',
    ha='left', va='center', fontsize=9, rotation=0
)
ax_amh.text(
    amh_primer_end + 600, -0.4,
    'Sysex-R: CCACCAATGATGGCGTCTAT',
    ha='right', va='center', fontsize=9, rotation=0
)

ax_amh.text(
    (amh_primer_start + amh_primer_end) / 2, 0.1,
    f'{amh_product_len} bp',
    ha='center', va='bottom', fontsize=12, fontweight='bold'
)

padding = (amh_end - amh_start) * 0.1
ax_amh.set_xlim(amh_start - padding, amh_end + padding)
ax_amh.set_ylim(-1, 1)
ax_amh.get_xaxis().set_major_formatter(mticker.FuncFormatter(lambda x, p: format(int(x), ',')))
ax_amh.set_yticks([])
ax_amh.spines['top'].set_visible(False)
ax_amh.spines['right'].set_visible(False)
ax_amh.spines['left'].set_visible(False)

# 成图 (Finalize Plot)
plt.tight_layout(pad=3.0)
plt.savefig("gene_map_primers.pdf")
plt.show()