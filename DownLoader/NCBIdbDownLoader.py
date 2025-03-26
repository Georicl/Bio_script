#!/bin/python
import os
import time
import hashlib
from urllib.parse import urljoin
import requests

BaseUrl = "https://ftp.ncbi.nlm.nih.gov/blast/db/"
FileNamePattern = "nt_prok.{0:2d}.tar.gz"
MD5NamePattern = "nt_prok.{0:2d}.tar.gz.md5"
OutputDir = "./"
MaxRetry = 1 # 重试次数


def Dowloadfile(Url, OutPutFile):
    "下载器"
    # 验证文件大小
    if os.path.exists(OutPutFile):
        CurrentSize = os.path.getsize(OutPutFile)
    else:
        CurrentSize = 0

    with requests.Session() as session:
        for attempt in range(MaxRetry):
            try:
                response = session.get(Url, timeout=30)
                if response.status_code == 416:
                    return True
                response.raise_for_status()

                # 检查文件大小
                Mode = "ab" if CurrentSize else "wb"

                with open(OutPutFile, Mode) as f:
                    # 每8MB写入一次块
                    for chunk in response.iter_content(chunk_size=8192):
                        if chunk:
                            f.write(chunk)
                return True
            except requests.exceptions.RequestException as e:
                print(f"Downloading happen an error in {OutPutFile}, retry ({attempt + 1}/{MaxRetry})")
                time.sleep(30) # 30s重试一次
        return False

def IdentifyMD5(FilePath, Md5Path):
    """Md5 校验器"""
    if not os.path.exists(Md5Path):
        return False

    with open(Md5Path, "r") as f:
        StdMD5 = f.read().strip().split()[0]

    MD5Hash = hashlib.md5()
    with open(FilePath, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            MD5Hash.update(chunk)

    FileMD5 = MD5Hash.hexdigest()
    return StdMD5 == FileMD5

def main():
    # 检验输出目录
    os.makedirs(OutputDir, exist_ok=True)

    # 获取文件列表
    FileList = [FileNamePattern.format(i) for i in range(0, 25)]
    MD5List = [MD5NamePattern.format(i) for i in range(0, 25)]

    # 下载文件
    for File, MD5 in zip(FileList, MD5List):
        FileUrl = urljoin(BaseUrl, File)
        MD5Url = urljoin(BaseUrl, MD5)
        FilePath = os.path.join(OutputDir, File)
        MD5Path = os.path.join(OutputDir, MD5)

        if IdentifyMD5(FilePath, MD5Path):
            continue

        Dowloadfile(FileUrl, FilePath)
        Dowloadfile(MD5Url, MD5Path)

    FailFile = []
    for File, MD5 in zip(FileList, MD5List):
        FilePath = os.path.join(OutputDir, File)
        MD5Path = os.path.join(OutputDir, MD5)

        if not IdentifyMD5(FilePath, MD5Path):
            FailFile.append(File)

if __name__ == "__main__":
    main()