# seqtools/__init__.py
import os
import importlib

__all__ = []
__path__ = [os.path.dirname(__file__)]

# 新增：排除工具脚本列表
ignore_files = ("_", "draw_plot", "cli")  # 新增"cli"过滤

for module_file in os.listdir(__path__[0]):
    if module_file.endswith(".py") and not module_file.startswith(ignore_files):
        module_name = module_file[:-3]
        module = importlib.import_module(f".{module_name}", package=__name__)

        for item in dir(module):
            if not item.startswith("_"):
                globals()[item] = getattr(module, item)
                __all__.append(item)
