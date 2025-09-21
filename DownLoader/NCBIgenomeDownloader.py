#!/usr/bin/env python3
import os
import time
import hashlib
import logging
import requests
from concurrent.futures import ThreadPoolExecutor
from urllib.parse import urljoin
import csv
import re
# 配置参数
BASE_URL = "https://ftp.ncbi.nlm.nih.gov/genomes/all/"
FILE_TYPES = {
    '_genomic.fna.gz': '.fna.gz',  # Genome assembly
    '_genomic.gff.gz': '.gff.gz',  # GFF annotations
    '_protein.faa.gz': '.faa.gz',  # Protein sequences
    '_cds_from_genomic.fna.gz': '.cds.fna.gz'  # Coding sequences
}
OUTPUT_DIR = "./downloads"
MAX_RETRY = 3
CHUNK_SIZE = 81920  # 下载缓存块大小，每80KB写入一次文件
THREADS = 5  # 并发下载数

# 初始化日志
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")


def normalize_filename(organism_name: str) -> str:
    """Sanitize organism names for filesystem safety"""
    return re.sub(r'[^\w-]', '_', organism_name).strip('_')


def construct_ftp_path(assembly_acc: str, assembly_name: str) -> str:
    """Build NCBI FTP path from assembly metadata"""
    try:
        prefix, identifier = assembly_acc.split('_', 1)
        base_acc = identifier.split('.', 1)[0]

        if len(base_acc) != 9:
            raise ValueError("Invalid assembly accession format")

        dir_structure = [
            prefix,
            base_acc[0:3],
            base_acc[3:6],
            base_acc[6:9],
            f"{assembly_acc}_{assembly_name}"
        ]

        return f"{BASE_URL}/{'/'.join(dir_structure)}/"

    except (ValueError, IndexError) as e:
        raise RuntimeError(f"Invalid assembly format: {assembly_acc}") from e


def download_file(url, output_file):
    """下载器，支持续传"""
    if os.path.exists(output_file):
        current_size = os.path.getsize(output_file)
    else:
        current_size = 0

    with requests.Session() as session:
        for attempt in range(MAX_RETRY):
            try:
                headers = {"Range": f"bytes={current_size}-"} if current_size else {}
                response = session.get(url, headers=headers, timeout=30, stream=True)
                response.raise_for_status()

                mode = "ab" if current_size else "wb"
                with open(output_file, mode) as f:
                    for chunk in response.iter_content(chunk_size=CHUNK_SIZE):
                        if chunk:
                            f.write(chunk)
                return True
            except requests.exceptions.RequestException as e:
                logging.error(f"Error downloading {output_file}: {e}, retry ({attempt + 1}/{MAX_RETRY})")
                time.sleep(30)  # 30秒后重试
        return False


def identify_md5(file_path, md5_path):
    """MD5 校验器"""
    if not os.path.exists(md5_path):
        logging.warning(f"MD5 file not found: {md5_path}")
        return False

    with open(md5_path, "r") as f:
        std_md5 = f.read().strip().split()[0]

    md5_hash = hashlib.md5()
    with open(file_path, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            md5_hash.update(chunk)

    file_md5 = md5_hash.hexdigest()
    return std_md5 == file_md5


def download_and_verify_assembly_files(row: list, file_types: dict):
    """Handle file downloading and verification logic for a single assembly"""
    assembly_acc = row[0].strip()
    assembly_name = row[1].strip()
    organism = row[2].strip()

    if not all([assembly_acc, assembly_name, organism]):
        logging.warning(f"Skipping incomplete record: {row}")
        return

    safe_name = normalize_filename(organism)
    logging.info(f"\nProcessing: {organism} ({assembly_acc})")

    try:
        ftp_base = construct_ftp_path(assembly_acc, assembly_name)
    except RuntimeError as e:
        logging.error(f"Construction error: {str(e)}")
        return

    for suffix, ext in file_types.items():
        remote_file = f"{assembly_acc}_{assembly_name}{suffix}"
        local_file = os.path.join(OUTPUT_DIR, f"{safe_name}{ext}")
        remote_url = urljoin(ftp_base, remote_file)

        if os.path.exists(local_file):
            logging.info(f"Existing file skipped: {local_file}")
            continue

        logging.info(f"Downloading {remote_url}")
        if download_file(remote_url, local_file):
            logging.info(f"Download complete: {local_file}")
        else:
            logging.error(f"Failed to download: {local_file}")
            continue

        md5_suffix = f"{suffix}.md5"
        md5_remote_file = f"{assembly_acc}_{assembly_name}{md5_suffix}"
        md5_local_file = os.path.join(OUTPUT_DIR, f"{safe_name}{md5_suffix}")
        md5_remote_url = urljoin(ftp_base, md5_remote_file)

        logging.info(f"Downloading MD5 checksum: {md5_remote_url}")
        if download_file(md5_remote_url, md5_local_file):
            logging.info(f"MD5 checksum downloaded: {md5_local_file}")
        else:
            logging.error(f"Failed to download MD5 checksum: {md5_local_file}")
            continue

        if identify_md5(local_file, md5_local_file):
            logging.info(f"MD5 verification successful: {local_file}")
        else:
            logging.error(f"MD5 verification failed: {local_file}")
            os.remove(local_file)
            os.remove(md5_local_file)


def main(input_file: str):
    """Main execution flow"""
    if not os.path.isfile(input_file):
        logging.error(f"Input file not found: {input_file}")
        sys.exit(1)

    os.makedirs(OUTPUT_DIR, exist_ok=True)

    with open(input_file, 'r') as tsv_file:
        tsv_reader = csv.reader(tsv_file, delimiter='\t')
        next(tsv_reader)  # Skip header

        with ThreadPoolExecutor(max_workers=THREADS) as executor:
            futures = []
            for record in tsv_reader:
                if len(record) >= 3:
                    futures.append(executor.submit(download_and_verify_assembly_files, record, FILE_TYPES))

            # 等待所有任务完成
            for future in futures:
                future.result()


if __name__ == '__main__':
    import sys

    if len(sys.argv) != 2:
        logging.error("Usage: python ncbi_assembly_downloader.py <input_file.tsv>")
        sys.exit(1)

    main(sys.argv[1])