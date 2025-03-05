# seqtools/__init__.py

# 显式导入所有子模块的函数
from .cli import seq_sort
from .cli import seq_shuffle
from .cli import seq_extract
from .cli import seq_bed_mapping
from .cli import subseq
from .cli import seq_search, read_seq
from .N50 import N50
from .Blast_identify import blast_identify
from .get_longest import get_longest
from .RNA_seq import RNA_seq
from .draw_plot import draw_plot
from .draw_plot import go_plot
from .draw_plot import run_plot

# 声明所有公共接口
__all__ = [
    'seq_sort',
    'seq_shuffle',
    'seq_extract',
    'seq_bed_mapping',
    'subseq',
    'seq_search',
    'read_seq',
    'RNA_seq',
    'N50',
    'blast_identify',
    'get_longest'
]
