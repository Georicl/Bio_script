#!/usr/bin/env python3
from collections import defaultdict


def parse_gff(gff_file):
    """
    解析GFF文件，返回按mRNA分组的CDS信息
    """
    mrna_dict = defaultdict(dict)

    with open(gff_file) as f:
        for line in f:
            if line.startswith('#') or line.strip() == '':
                continue

            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue

            chrom, source, feat_type, start, end, score, strand, phase, attributes = fields

            if feat_type not in ('mRNA', 'CDS'):
                continue

            attr_dict = {}
            for attr in attributes.split(';'):
                if '=' in attr:
                    key, value = attr.split('=', 1)
                    attr_dict[key.strip()] = value.strip()

            if feat_type == 'mRNA':
                mrna_id = attr_dict.get('ID')
                if mrna_id:
                    mrna_dict[mrna_id] = {
                        'start': int(start),
                        'end': int(end),
                        'cds': []
                    }
            elif feat_type == 'CDS':
                parent_id = attr_dict.get('Parent')
                if parent_id and parent_id in mrna_dict:
                    mrna_dict[parent_id]['cds'].append(
                        (int(start), int(end))
                    )

    return mrna_dict


def calculate_top_introns(mrna_data, top_n=3):
    """
    计算内含子长度并返回前N名
    """
    intron_lengths = []

    for mrna_id, data in mrna_data.items():
        mrna_length = data['end'] - data['start'] + 1
        cds_total = sum(end - start + 1 for (start, end) in data['cds'])
        intron_length = mrna_length - cds_total

        # 只保留有效非负长度
        if intron_length >= 0:
            intron_lengths.append(intron_length)

    # 降序排序并取前N项
    sorted_lengths = sorted(intron_lengths, reverse=True)
    return sorted_lengths[:top_n]


if __name__ == '__main__':
    import sys

    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} <input.gff>")
        sys.exit(1)

    mrna_data = parse_gff(sys.argv[1])
    top_introns = calculate_top_introns(mrna_data)

    print("Top 3 intron lengths:")
    for i, length in enumerate(top_introns, 1):
        print(f"{i}. {length} bp")

    # 补充输出有效内含子总数
    print(f"\nTotal valid introns: {len(top_introns)}")
