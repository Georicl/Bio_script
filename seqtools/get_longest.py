import os
import gffutils
from Bio import SeqIO


def get_longest(cds,out_path,gff):
    #########################
    #the function for get the longest transcripts in cds sequence
    #input cds sequence  out_path gff3 file path
    #out_path : /path/to/file.fasta
    #return none
    ########################
    cdsseq = SeqIO.index(cds,'fasta')
    gffdb = gffutils.create_db(gff,'gff.db',force=True,merge_strategy='create_unique')
    gene_to_longestcds = {}
    out = os.path.jion(out_path,'longest_cds.fasta')
    for g in gffdb.all_features(featuretype='gene'):
        g_id = g.id
        for m in gffdb.children(g,featuretype='mRNA'):
            m_id = m.id.replace('transcript:', '')
            if m_id not in cdsseq:
                print('transcript'+ m_id + 'not in cds')
                continue
            else:
                m_len = len(cdsseq[m_id].seq)
            if g_id not in gene_to_longestcds:
                gene_to_longestcds[g_id] = m_id
            elif m_len > len(cdsseq[gene_to_longestcds[g_id]].seq):
                gene_to_longestcds[g_id] = m_id
    sr_list = [v for k,v in cdsseq.items() if k in gene_to_longestcds.values()]
    SeqIO.write(sr_list,out,'fasta')
    return