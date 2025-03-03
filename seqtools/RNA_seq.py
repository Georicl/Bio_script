import subprocess

def RNA_seq(RNA_Seq_path, fasta_path, gtf_path, sample_file, cpu, method, data_path, script_path, work_path):
    cmd = f"sh {RNA_Seq_path} {fasta_path} {gtf_path} {sample_file} {cpu} {method} {data_path} {script_path} {work_path}"
    subprocess.run(cmd, shell=True, check=True)