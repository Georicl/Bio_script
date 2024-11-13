#######################################
#region import package
######################################
import argparse
import math
import random
import os
import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import Blast_identify
import N50
import get_longest


#endregion
######################################
###region function
######################################
###random gene extract
def seq_extract(fasta_path,number,out_path,seed=10):
    ################################################
    #function for extract genes random
    ################################################
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
    return

def seq_sort(fasta_path,out_path,reverse):
    ################################################
    #function for sort fasta
    ################################################
    out_sort_path = os.path.join(out_path, 'sort.fasta')  # get out path
    seq_index = SeqIO.index(fasta_path,'fasta') #index gene fasta
    seq_sort_list = sorted(list(seq_index.items()),key=lambda x:len(x[1].seq),reverse=reverse) #sort gene by genes' length
    ##output
    SeqIO.write(seq_sort_list,out_sort_path,'fasta')
    return

def seq_shuffle(fasta_path,out_path):
    ####################################
    #shuffle function
    ####################################
    out_of_order = os.path.join(out_path, 'shuffle.fasta')  # get out path
    seq_index = SeqIO.index(fasta_path, 'fasta')  # index gene fasta
    out_of_order_list = list(seq_index.values())
    random.shuffle(out_of_order_list) #shuffle
    #output
    SeqIO.write(out_of_order_list,out_of_order,'fasta')
    return

def subseq(fasta_path,sub_data,out_path):
    ################################################
    #sub_seq by id
    ################################################
    out_seq_path = os.path.join(out_path, 'subseq.fasta')  # get out path
    seq_index = SeqIO.index(fasta_path, 'fasta')  # index gene fasta
    with open(sub_data,'r') as finput: #read id list
        read_keys = finput.read().strip().split('\n')
        read_list = [v for k,v in seq_index.items() if k in read_keys ] #mapping
    #output
    SeqIO.write(read_list,out_seq_path,'fasta')
    return

def seq_bed_mapping(fasta_path,bed_data,out_path):
    ##########################################
    #sub_seq by bed
    ##########################################
    out_bed_path = os.path.join(out_path, 'bed.fasta')  # get out path
    seq_index = SeqIO.index(fasta_path, 'fasta')  # index gene fasta
    seq_record_list = []
    with open(bed_data,'r') as finput: #read bed
        for line in finput:
            parts = line.strip().split("\t")
            acc = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            #extract all genes cds
            sequence = seq_index[acc].seq
            #extract the sequence you want to use
            sub_sequence = sequence[start-1:end]
            #build SeqRecord
            seq_record = SeqRecord(id=':'.jion((acc,str(start)),str(end)),
                                        seq=sub_sequence,
                                   description='seq_index[acc].description')
            seq_record_list.append(seq_record)
            #output
            SeqIO.write(seq_record_list,out_bed_path,'fasta')
            return
#############################
def read_seq(search_input):
    ########################
    #read sequence
    #return the sequence
    ########################
    with open(search_input,'r') as finput:
        sequence = finput.read()
    return sequence

def seq_search(fasta_path,search_input,out_path):
    ##########################################
    #sub_fasta by sequence
    ##########################################
    out_sequence_path = os.path.join(out_path, 'searchid.fasta')  # get out path
    sr_list = []
    seq_index = SeqIO.index(fasta_path, 'fasta')  # index gene fasta
    for k,v in seq_index.items():
        sequence = str(v.seq)
        start = sequence.find(search_input)
        SeqRecord(id=seq_index[k].id,seq=seq_index[k].seq,description=seq_index[k].description+'_'+str(start))
    SeqIO.write(sr_list,out_sequence_path,'fasta')
    return

##############################
#get hisat2.log
##############################
def read_hisat2_log(sample_txt, input_path,out_path):
    ########################################
    #function for filter hisat2.log
    #only extracted contains alignment rate percent
    ########################################
    output_ht2log_path = os.path.join(out_path, 'ht2log.txt')  # get out path
    with open(sample_txt, 'r') as sample_file:
        filenames_temp = sample_file.readlines()
        filenames = []
        for line in filenames_temp: #read sample_txt and extract the first colunmn
            line_data = line.strip().split('\t')
            filenames.append(line_data[0])
        results = {}
        line_need = []
        for line in filenames:
            full_path = os.path.join(input_path, line + '.log')# get full input path
            line_name = line
            with open(full_path) as f:
                for line_1 in f:
                    line_1 = line_1.strip()
                    if line_1.startswith("Overall"): #filter rows which start with "Overall"
                        line_need = line_1.split(':')[1]
                        #line_need is percent of alignment rate in the log
            results[line_name] = line_need
    with open(output_ht2log_path, 'w') as output_file:
        for key, value in results.items():
            output_file.write(f"{key}\t{value}\n")
    return

#############################
#endregion
######################################
######################################
###main
######################################
###region parse argument
#create arg

parse = argparse.ArgumentParser(description = "This tool is used to sequence processing, fasta"
                                              " reading and information extraction, and can also be used to process the "
                                              "log files generated by blast and hisat2 and for quality checking of log files,"
                                              "and calculate N50 and GC content ")
subparsers = parse.add_subparsers(dest='command',required=True,help='subcommands')

sort_parser = subparsers.add_parser('sort', help='Sort sequences',usage='seqtools.py sort [options] gene_fasta out_path')
sort_parser.add_argument('gene_fasta', type=str, help='Input your genome FASTA file')
sort_parser.add_argument('out_path', type=str, help='Output path directory')
sort_parser.add_argument('-r','--reverse', action='store_true', help='Reverse sort order',default=False)

shuffle_parser = subparsers.add_parser('shuffle',help = 'Shuffle sequence')
shuffle_parser.add_argument('gene_fasta', type=str, help='Input your genome FASTA file')
shuffle_parser.add_argument('out_path', type=str, help='Output path directory')

random_parser = subparsers.add_parser('extract',help='Extract random sequences')
random_parser.add_argument('gene_fasta', type=str, help='Input your genome FASTA file')
random_parser.add_argument('out_path', type=str, help='Output path directory')
random_parser.add_argument('number',type=float,help='input your number of gene you want to extract')

bed_parse = subparsers.add_parser('bed',help='Filter BED file')
bed_parse.add_argument('gene_fasta', type=str, help='Input your genome FASTA file')
bed_parse.add_argument('out_path', type=str, help='Output path directory')
bed_parse.add_argument('bed',type=str,help='input your bed file')

subseq_parse = subparsers.add_parser('subseq',help='sub sequences by geneid')
subseq_parse.add_argument('gene_fasta', type=str, help='Input your genome FASTA file')
subseq_parse.add_argument('out_path', type=str, help='Output path directory')
subseq_parse.add_argument('sub_data',type=str,help='input your sub data file')

search_parse = subparsers.add_parser('search',help='Search subsequences by Nucleotide sequences')
search_parse.add_argument('gene_fasta', type=str, help='Input your genome FASTA file')
search_parse.add_argument('out_path', type=str, help='Output path directory')
search_parse.add_argument('sequence',type=str,
                          help='input sequence which you want to extract')

n50_parse = subparsers.add_parser('N50',help='calculate input fasta file and get N50 and GC content，'
                                             'Statistical series parameters')
n50_parse.add_argument('gene_fasta',type=str, help='Input your genome FASTA file')
n50_parse.add_argument('out_path', type=str, help='Output path directory')

blast_parse = subparsers.add_parser('blast',help='filter indentify you need,if you need to use this mod,please use '
                                                 'blast outputformat 6')
blast_parse.add_argument('blast_in', type=str, help='Input your blast file')
blast_parse.add_argument('out_path', type=str, help='Output path directory')
blast_parse.add_argument('indentify', type=float, help='Filter how much recognition rate is above ,'
                                                       'The input units are %')

get_longest_parse = subparsers.add_parser('get_longest',help='get longest cds sequence')
get_longest_parse.add_argument('gene_fasta', type=str, help='Input your genome FASTA file')
get_longest_parse.add_argument('out_path', type=str, help='Output path directory')
get_longest_parse.add_argument('gff',type=str,help='input gff3 file')
#endregion
##################################################
#region for hisat2.log
##################################################
ht2log_parse = subparsers.add_parser('ht2log',help='Process hisat2.log file')
ht2log_parse.add_argument('gene_fasta', type=str, help='Input your genome FASTA file')
ht2log_parse.add_argument('out_path', type=str, help='Output path directory')
ht2log_parse.add_argument('ht2log',type=str,help='input your hisat2.log directory,don\'t have / in tail')
ht2log_parse.add_argument('sample',type=str,help='input your sample file')
#endregion
######################################
##region help out
######################################
if len(sys.argv) == 1:
    parse.print_help()
    sys.exit(1)
if len(sys.argv) == 2 and sys.argv[1] =='sort':
    sort_parser.print_help()
    sys.exit(1)
elif len(sys.argv) == 2 and sys.argv[1] =='shuffle':
    shuffle_parser.print_help()
    sys.exit(1)
elif len(sys.argv) == 2 and sys.argv[1] =='extract':
    random_parser.print_help()
    sys.exit(1)
elif len(sys.argv) == 2 and sys.argv[1] =='bed':
    bed_parse.print_help()
    sys.exit(1)
elif len(sys.argv) == 2 and sys.argv[1] =='subseq':
    subseq_parse.print_help()
    sys.exit(1)
elif len(sys.argv) == 2 and sys.argv[1] =='search':
    search_parse.print_help()
    sys.exit(1)
elif len(sys.argv) == 2 and sys.argv[1] =='ht2log':
    ht2log_parse.print_help()
    sys.exit(1)
elif len(sys.argv) == 2 and sys.argv[1] =='N50':
    n50_parse.print_help()
    sys.exit(1)
elif len(sys.argv) == 2 and sys.argv[1] =='blast':
    blast_parse.print_help()
    sys.exit(1)
#endregion
######################################
#parse
######################################
args = parse.parse_args()
######################################
#main function
######################################
if args.command =='sort':
    seq_sort(args.gene_fasta,args.out_path,args.reverse)
elif args.command =='shuffle':
    seq_shuffle(args.gene_fasta,args.out_path)
elif args.command =='extract':
    seq_extract(args.gene_fasta,args.number,args.out_path)
elif args.command =='bed':
    seq_bed_mapping(args.gene_fasta,args.bed_data,args.out_path)
elif args.command =='subseq':
    subseq(args.gene_fasta,args.sub_data,args.out_path)
elif args.command =='search':
    sequence = read_seq(args.search_input)
    seq_search(args.gene_fasta,sequence,args.out_path)
elif args.command =='ht2log':
    read_hisat2_log(args.sample, args.ht2log,args.out_path)
############################
#import module
############################
elif args.command == 'N50': #import N50
    N50.N50(args.gene_fasta,args.out_path)
elif args.command == 'blast':
    Blast_identify.blast_identify(args.blast_in, args.out_path, args.indentify)
elif args.command == 'get_longest':
    get_longest.get_longest(args.cds, args.out_path, args.gff)
#end