import argparse
import os 
import requests
import configparser
import hashlib
from concurrent.futures import ThreadPoolExecutor, as_completed  # 并发模块

class DownLoader:
    """
    下载器，使用配置文件config.ini实现对不同url，以及不同url中的需要下载的文件解析。
    支持使用tsv中第一列为下载清单列表。
    支持使用直接输入的sample = sample1, sample2的方法输入多个文件
    支持输入多个url如使用：[url]的方式间隔区分。
    支持多重下载。
    """
    def __init__(self, config_path):
        # 修改参数
        self.retry = 5 # 重试次数
        self.chunk_size = 81920 # 缓存大小80kb
        self.threads = 5 # 并发下载数
        self.output_path = "./"
        # 读取配置文件
        self.config_file = self._parser_config(config_path)

    def _tsv_process_sample(self, sample_tsv_path: str, http_url: str):
        sample_list = []
        with open(sample_tsv_path, 'r') as f:
            for lines in f:
                filename = lines.strip().split('\t')[0]
                sample_list.append({'url': f"{http_url}/{filename}",
                                    'filename': filename,
                                    'md5_url': f"{http_url}/{filename}.md5"
                                    })
        return sample_list
    
    def _process_sample(self, sample, http_url: str):
        # 判断tsv
        if isinstance(sample, str) and sample.endswith(".tsv"):
            return self._tsv_process_sample(sample, http_url)
        
        if isinstance(sample, str):
            sample = [s.strip().strip("'\"") for s in sample.split(',')]

        sample_list = []
        for filename in sample:
            sample_list.append({'url': f"{http_url}/{filename}",
                                'filename': filename,
                                'md5_url': f"{http_url}/{filename}.md5"
                                })
        return sample_list
    
    def _parser_config(self, config_path):
        # 解析config
        config = configparser.ConfigParser()
        config.read(config_path)
        
        # 解析下载请求的路径的文件
        download_list = []
        for section in config.sections():
            if section.startswith('url'):
                config_file = {
                    'http_url': config[section]['url'].strip('""'),
                    'md5_value': config[section].getboolean('md5', False),
                    'samples': self._process_sample(config[section]['sample'], 
                                                   config[section]['url'].strip('"')
                                                   )
                }
            download_list.append(config_file)
        return download_list

    def _retry_download(self, url, download_path):
        for times in range(self.retry + 1):
            try:
                # 请求下载
                with requests.get(url, stream=True) as r:
                    r.raise_for_status()
                    with open(download_path, 'wb') as f:
                        for chunk in r.iter_content(chunk_size=self.chunk_size):
                            f.write(chunk)
            except:
                if times == self.retry:
                    raise RuntimeError(f"Failed to download {url} after {self.retry} attempts")
                continue

    def _verify_md5(self, file_path, md5_url):
        # md5校验器
        expected_md5 = requests.get(md5_url).text.split()[0]
        hash_md5 = hashlib.md5()
        with open(file_path, "rb") as f:
            for chunk in iter(lambda: f.read(4096), b""):
                hash_md5.update(chunk)
        actual_md5 = hash_md5.hexdigest()
        if actual_md5 != expected_md5:
            raise ValueError(f"MD5 mismatch for {file_path}")
    
    def _download_task(self, sample, task):
        # 单个下载模块
        local_path = os.path.join(self.output_path, sample['filename'])
        try:
            self._retry_download(sample['url'], local_path)
            
            if task['md5_value']:
                md5_path = local_path + '.md5'
                self._retry_download(sample['md5_url'], md5_path)
                self._verify_md5(local_path, sample['md5_url'])
                
            with self.lock:
                print(f"成功下载 {sample['filename']}")
            return True
        except Exception as e:
            with self.lock:
                print(f"下载失败 {sample['filename']};原因: {str(e)}")
            return False

    def main(self):
        # 线程池
        with ThreadPoolExecutor(max_workers=self.threads) as executor:
            futures = []
            for task in self.config_file:
                os.makedirs(self.output_path, exist_ok=True)
                
                # 提交并发任务
                for sample in task['samples']:
                    future = executor.submit(
                        self._download_task, 
                        sample=sample,
                        task=task
                    )
                    futures.append(future)
                    
            # 打印任务成功数量
            success = 0
            failure = 0
            for future in as_completed(futures):
                if future.result():
                    success += 1
                else:
                    failure += 1
            
            print(f"\n下载完成 | 成功: {success} 个 | 失败: {failure} 个")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="下载器，根据指定url批量下载样本")
    parser.add_argument("config_path", type=str, help="指定配置文件")
    args = parser.parse_args()
    # 初始化下载管理器并执行任务
    downloader = DownLoader(args.config_path)
    downloader.main()
