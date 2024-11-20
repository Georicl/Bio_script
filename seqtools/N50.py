import os
import sys
from pathlib import Path
import re
import pysam
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
##################################
#open fa file used seq
def N50(fa_argv,out_path,trinity_gene):
    ###################################
    #function for calculate N50 and total length
    ###################################
    ###################################
    #create fai file
    ###################################
    pysam.faidx(fa_argv)
    fai = fa_argv+'.fai'
    seq_index= SeqIO.index(fa_argv,'fasta')
    #chr_len = {k:len(v.seq) for k,v in seq_index.items()}
    out_N50_path = os.path.join(out_path,'N50.txt')
    #open fai file
    with open(fai, 'r') as f:
        chr_len = {}
        for line in f:
            info = line.strip().split('\t')
            chr_len[info[0]] = int(info[1])
    ###################################
    ### endregion
    ###################################
    ### region calculate
    #n_count,gaps
    gaps = 0
    n_count = 0
    total_gc = 0
    for k,v in seq_index.items():
        s = str(v.seq).upper() #convert to upper
        n_count += s.count('N')
        gaps += len(re.split('N+',s))-1
        total_gc += s.count('G') + s.count('C') #calculate GC

    #sum length ,n, mean_n
    n = len(chr_len)
    total_n = sum(chr_len.values())
    mean_n = total_n / n
    GC_content = total_gc / total_n * 100 if total_n > 0 else 0
    ###################################
    ###################################
    #longest sequence
    #longest_id, longest_chr = max(chr_len.items(),key = lambda  x: x[1])

    ###################################
    def nx_sum(nx_long, nx):
        #############################################
        ###function for Calculate the length of n(x).
        ###use nx_long,nx,
        ###return nx1:total length for n(x)
        ###       idx:the number of that sequence
        #############################################
        #sort and Initialization parameters
        nx_long = sorted(nx_long, reverse=True)
        nx1 = 0
        idx = 0
        sum_nx = 0
        p = int(nx[1:]) / 100
        for i,v in enumerate(nx_long, start = 1):
            sum_nx += v
            if sum_nx >= total_n*p :
                nx1 = v
                idx = i
                break
        return nx1,idx
    ##################################################
    # endregion
    ##################################################
    path_fa = Path(fa_argv)
    if trinity_gene == 1:
        ####################
        #if input is a trinity fasta and want to statistic
        ####################
        with open(fai,'r') as f:
            numgenes = 0
            for line in f:
                info = line.strip().split('\t')[0]
                pattern = r'_\d+_'
                matchs = re.findall(pattern,info)
                matches = re.findall(r'\d+',matchs[0])
                if int(matches[0]) > numgenes:
                    numgenes = int(matches[0])
    with open(out_N50_path,'w') as f:
    #write results
        f.write(f'File={path_fa.stem}\n'
                f'Total length = {total_n:<10}\n'
                f'Mean length = {mean_n:<10}\n'
                f'GC content = {GC_content:<10}\n')  #only print file name
        for n in ('N50','N60','N70','N80','N90','N100'): #print N50,N60,N70,N80,N90,N100
                nx1,idx = nx_sum(chr_len.values(), n)
                f.write(f'{n:<4} = {nx1:<10}, n = {idx:<10}\n')
        f.write(f'N count= {n_count:<10}, Gaps = {gaps:<10}\n')
        if trinity_gene == 1:
            ####################
            # if input is a trinity fasta and want to statistic
            ####################
            f.write(f'This is trinity gene argument:\n'
                f'  gene = {numgenes:<10}\n')

    return