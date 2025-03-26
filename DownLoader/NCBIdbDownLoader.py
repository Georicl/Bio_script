#!/usr/bin/env python3
import os
import time
import hashlib
import logging
from urllib.parse import urljoin
import requests
from concurrent.futures import ThreadPoolExecutor

# 配置参数
BASE_URL = "https://ftp.ncbi.nlm.nih.gov/blast/db/"
FILE_NAME_PATTERN = "nt_prok.{}.tar.gz"
MD5_NAME_PATTERN = "nt_prok.{}.tar.gz.md5"
OUTPUT_DIR = "./"
MAX_RETRY = 3
CHUNK_SIZE = 8192  # 下载缓存块大小
THREADS = 5  # 并发下载数

# 初始化日志
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

def download_file(url, output_file):
    """下载器"""
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

def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    file_list = [FILE_NAME_PATTERN.format(i) for i in range(25)]
    md5_list = [MD5_NAME_PATTERN.format(i) for i in range(25)]

    fail_files = []

    with ThreadPoolExecutor(max_workers=THREADS) as executor:
        futures = []
        for file, md5 in zip(file_list, md5_list):
            file_url = urljoin(BASE_URL, file)
            md5_url = urljoin(BASE_URL, md5)
            file_path = os.path.join(OUTPUT_DIR, file)
            md5_path = os.path.join(OUTPUT_DIR, md5)

            futures.append(executor.submit(download_file, md5_url, md5_path))
            futures.append(executor.submit(download_file, file_url, file_path))

        # 等待所有任务完成
        for future in futures:
            future.result()

    # 校验文件
    for file, md5 in zip(file_list, md5_list):
        file_path = os.path.join(OUTPUT_DIR, file)
        md5_path = os.path.join(OUTPUT_DIR, md5)

        if not identify_md5(file_path, md5_path):
            fail_files.append(file)

    if fail_files:
        logging.error(f"Failed to verify the following files: {fail_files}")
    else:
        logging.info("All files downloaded and verified successfully.")

if __name__ == "__main__":
    main()
