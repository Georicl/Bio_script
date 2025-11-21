#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pickle
import csv
import argparse


def convert_lefse_to_csv(lefse_in_path, sample_map_path, output_csv_path):
    """
    将 LEfSe 的 .in 输入文件和样本映射文件转换为标准的 CSV 表格。

    Args:
        lefse_in_path (str): LEfSe .in 文件的路径 (pickle 格式)。
        sample_map_path (str): 样本映射文件的路径 (第一列样本ID，第二列分组ID)。
        output_csv_path (str): 输出 CSV 文件的路径。
    """
    print("🚀 开始转换 LEfSe 文件...")

    # --- 1. 读取样本和分组信息 ---
    # 这个顺序决定了数据列的顺序
    sample_ids = []
    group_ids = []
    try:
        with open(sample_map_path, 'r') as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 2:
                    sample_ids.append(parts[0])
                    group_ids.append(parts[1])
        print(f"✅ 成功读取 {len(sample_ids)} 个样本信息。")
    except FileNotFoundError:
        print(f"❌ 错误: 找不到样本映射文件 '{sample_map_path}'")
        return

    # --- 2. 加载 lefse.in 文件 ---
    try:
        with open(lefse_in_path, 'rb') as f:
            # 使用 pickle.load() 来解析二进制文件
            lefse_data = pickle.load(f)
        print("✅ 成功加载 lefse.in 文件。")
    except FileNotFoundError:
        print(f"❌ 错误: 找不到 LEfSe 输入文件 '{lefse_in_path}'")
        return
    except pickle.UnpicklingError:
        print(f"❌ 错误: 文件 '{lefse_in_path}' 不是一个有效的 pickle 文件。")
        return

    # 从加载的数据中提取物种丰度信息
    features_data = lefse_data.get('feats', {})
    if not features_data:
        print("❌ 错误: lefse.in 文件中未找到 'feats' 数据。")
        return

    # --- 3. 写入 CSV 文件 ---
    with open(output_csv_path, 'w', newline='', encoding='utf-8') as csvfile:
        writer = csv.writer(csvfile)

        # 写入第一行表头：分组信息
        # 第一个单元格留空或写 'Group'，后面是每个样本对应的分组
        writer.writerow(['Group'] + group_ids)

        # 写入第二行表头：样本ID
        # 第一个单元格写 'Taxon' 或 'Feature'，后面是具体的样本ID
        writer.writerow(['Feature'] + sample_ids)

        # 按字母顺序写入每个物种（Feature）及其丰度数据
        # 排序可以确保每次运行输出的行顺序一致
        sorted_features = sorted(features_data.keys())

        for feature_name in sorted_features:
            abundances = features_data[feature_name]
            # 确保丰度列表长度与样本数量匹配
            if len(abundances) == len(sample_ids):
                # 将物种名称和其对应的丰度列表写入一行
                writer.writerow([feature_name] + abundances)
            else:
                print(
                    f"⚠️ 警告: 物种 '{feature_name}' 的数据点数量 ({len(abundances)}) 与样本数量 ({len(sample_ids)}) 不匹配，已跳过。")

    print(f"🎉 转换完成！输出文件已保存至: {output_csv_path}")


if __name__ == '__main__':
    # --- 设置命令行参数解析 ---
    parser = argparse.ArgumentParser(
        description="将 LEfSe 输入文件 (.in) 转换为人类可读的 CSV 格式。",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument('-i', '--input', required=True, help="输入的 lefse.in 文件路径。")
    parser.add_argument('-s', '--samples', required=True,
                        help="样本信息文件路径 (Tab或空格分隔, 第1列样本ID, 第2列分组ID)。")
    parser.add_argument('-o', '--output', required=True, help="输出的 CSV 文件路径。")

    args = parser.parse_args()

    # --- 调用主函数 ---
    convert_lefse_to_csv(args.input, args.samples, args.output)