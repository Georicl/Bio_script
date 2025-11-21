import pickle
import pandas as pd


def convert_lefse_to_tsv(lefse_file_path, sample_map_path, output_file_path):
    """
    将 LEfSe 的输入文件 (.in) 转换为标准的 TSV 表格文件。

    参数:
    lefse_file_path (str): LEfSe 输入文件 (.in) 的路径。
    sample_map_path (str): 样本映射文件的路径 (第一列: 样本ID, 第二列: 分组)。
    output_file_path (str): 输出的 TSV 文件路径。
    """
    print("--- 开始转换 ---")

    # 读取样本信息文件，获取有序的样本ID列表
    try:
        with open(sample_map_path, 'r') as f:
            sample_ids = [line.strip().split()[0] for line in f if line.strip()]
        print(f"成功读取 {len(sample_ids)} 个样本ID，从: {sample_map_path}")
    except FileNotFoundError:
        print(f"错误: 找不到样本映射文件 '{sample_map_path}'")
        return
    except Exception as e:
        print(f"读取样本映射文件时出错: {e}")
        return

    # 读取 LEfSe 的 pickle 输入文件
    try:
        with open(lefse_file_path, 'rb') as f:
            data = pickle.load(f)
        print(f"成功解析 LEfSe 输入文件: {lefse_file_path}")
    except FileNotFoundError:
        print(f"错误: 找不到 LEfSe 输入文件 '{lefse_file_path}'")
        return
    except pickle.UnpicklingError:
        print(f"错误: 文件 '{lefse_file_path}' 不是一个有效的 pickle 文件。")
        return
    except Exception as e:
        print(f"读取 LEfSe 文件时出错: {e}")
        return

    # 提取特征数据
    # 'feats' 是字典，其中键是物种/功能名，值是丰度列表
    if 'feats' in data and isinstance(data['feats'], dict):
        feature_data = data['feats']
        print(f"文件包含 {len(feature_data)} 个特征。")
    else:
        print("错误: LEfSe 文件中未找到 'feats' 键或其格式不正确。")
        return

    # 验证样本数量是否匹配
    # 随机抽取一个特征来检查其丰度列表的长度
    first_feature_name = next(iter(feature_data))
    num_values = len(feature_data[first_feature_name])
    if len(sample_ids) != num_values:
        print(f"警告: 样本ID数量 ({len(sample_ids)}) 与数据点数量 ({num_values}) 不匹配！")
        print("请检查你的样本映射文件是否与生成 lefse.in 文件时所用的数据一致。")
        # 如果列名不匹配, 可能产生列名错误，但可以继续运行

    # 创建 DataFrame
    # orient='index' 表示字典的键作为行索引
    print("正在创建数据表...")
    df = pd.DataFrame.from_dict(feature_data, orient='index')

    # 将列名设置为我们从样本文件中读取的ID
    df.columns = sample_ids

    # 保存Dataframe 为 TSV
    try:
        df.to_csv(output_file_path, sep='\t', index_label='Feature')
        print(f"--- 转换成功 ---")
        print(f"标准表格文件已保存至: {output_file_path}")
    except Exception as e:
        print(f"保存文件时出错: {e}")


if __name__ == "__main__":
    # 定义输入和输出文件路径
    LEFSE_INPUT_FILE = 'lefse.in'
    SAMPLE_MAPPING_FILE = 'sample_mapping.txt'  # 请确保你的样本信息文件名正确
    OUTPUT_TSV_FILE = 'lefse_output.tsv'

    # 调用转换函数
    convert_lefse_to_tsv(LEFSE_INPUT_FILE, SAMPLE_MAPPING_FILE, OUTPUT_TSV_FILE)