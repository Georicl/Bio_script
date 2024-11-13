import math
import os
from random import random
import argparse
from Bio import SeqIO


def seq_extract(fasta_path, number, out_path, seed=10,):
    ################################################
    #function for extract genes random
    ################################################
    if args.seed:
        seed = int(args.seed)
        out_extract_path = os.path.join(out_path, 'extract.fasta')
        seq_index = SeqIO.index(fasta_path,'fasta') #index gene fasta
        if number < 1 :
            number = math.ceil(len(seq_index)*number) #get number use percent
        elif number >= 1 :
            number = int(number)  #if number > 1,get int number
        random.seed(seed) #set seed
        keys = random.sample(list(seq_index.keys()),number) #random gene
        seq_index_out = {k:v for k,v in seq_index.items() if k in keys} #get random gene
        #ouput
        SeqIO.write(seq_index_out.values(),out_extract_path,'fasta')
    else:
        seed = int(args.seed)
        out_extract_path = os.path.join(out_path, 'extract.fasta')
        seq_index = SeqIO.index(fasta_path, 'fasta')  # index gene fasta
        if number < 1:
            number = math.ceil(len(seq_index) * number)  # get number use percent
        elif number >= 1:
            number = int(number)  # if number > 1,get int number
        random.seed(seed)  # set seed
        keys = random.sample(list(seq_index.keys()), number)  # random gene
        seq_index_out = {k: v for k, v in seq_index.items() if k in keys}  # get random gene
        # ouput
        SeqIO.write(seq_index_out.values(), out_extract_path, 'fasta')
    return