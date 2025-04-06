from os import path


def get_paths( f_bam, m_bam, fasta_path_get, output_path_get_path):
    """
    Generate and return dictionary of output paths.
    Input: fasta
    Output: file path
    """
    temp_dir = path.join(output_path_get_path, "temp_dir")
    return {
        "fa_path": f"{fasta_path_get}",
        "log_path": path.join(output_path_get_path, "CQ.log"),
        "fai_path": f"{fasta_path_get}.fai",
        "chromosome_length": path.join(output_path_get_path, "chromosome_length.tsv"),
        "windows_tsv": path.join(output_path_get_path, "chromosome_windows.tsv"),
        "f_cov": path.join(temp_dir, "F_windows_coverage.tsv"),
        "m_cov": path.join(temp_dir, "M_windows_coverage.tsv"),
        "f_m_merge": path.join(temp_dir, "F_M_merge.tsv"),
        "merged_result": path.join(temp_dir, "F_M_merge_result.tsv"),
        "cq_output_f": path.join(temp_dir, "F_M_CQ.tsv"),
        "reads_file": path.join(temp_dir, "F_M_reads"),
        "filtered_cq": path.join(output_path_get_path, "F_M_CQ.filter.tsv"),
        "tmm_normalized": path.join(temp_dir, "F_M_reads_TMM.tsv"),
        "temp_dir": path.join(output_path_get_path, "temp_dir"),
        "f_bam": f_bam,
        "m_bam": m_bam,
    }