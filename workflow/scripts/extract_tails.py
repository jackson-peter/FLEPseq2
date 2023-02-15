#!/usr/bin/env python3
import pandas as pd
import numpy as np
import swifter
import click
import pyfastx
import sys
from Bio import SeqIO
import gzip
import edlib
from itertools import zip_longest
import regex
from tqdm import tqdm
from colorama import Fore, Style, init
from collections import Counter
init(autoreset=True)


equalities=[("M", "A"), ("M", "C"),("R", "A"), ("R", "A"), ("W", "A"), ("W", "A"), ("S", "C"), ("S", "C"), ("Y", "C"), ("Y", "C"), 
("K", "G"), ("K", "G"), ("V", "A"), ("V", "C"), ("V", "G"), ("H", "A"), ("H", "C"), ("H", "T"), ("D", "A"), ("D", "G"), ("D", "T"),
("B", "C"), ("B", "G"), ("B", "T"), ("N", "G"), ("N", "A"), ("N", "T"), ("N", "C")]
    

FWD=("F-R", "F-UR", "F-UUR", "F-UF", "F-F", "F-N", "N-R", "N-N", "N-UR", "N-UUR", "UF-R", "UF-N", "UF-UF", "UF-UR", "UF-UUR", "UR-R", "UUR-R", "UUR-UR")
REV=("R-F", "R-R", "R-N", "R-UF", "R-UUR", "R-UR", "N-F", "N-UF", "UF-F", "UR-F", "UR-N", "UR-UF", "UR-UR", "UR-UUR", "UUR-F", "UUR-N", "UUR-UF", "UUR-UUR")
 
@click.command()
@click.option('-i', '--inadapter', help='Input adapter file', required=True, type=click.Path(exists=True))                         
@click.option('-s', '--inseq', help='Input read fastq file', required=True, type=click.Path(exists=True))                             
@click.option('-o', '--out', help='Output adapter information of each read', required=True)
@click.option('-v', '--verbose', is_flag=True, help="Print sequences with colors for polytail, adapter and delimiter")
@click.option('-d', '--debug', is_flag=True, help="for developping purposes, prints additional information")

def main(inadapter, inseq, out, constant_seq="CTGAC", umi_seq="NNNNNNNNNN", adapt_seq="CTGTAGGCACCATCAAT", verbose=False, debug=False):
    fields = ['read_core_id','mRNA','init_polya_start_base', 'init_polya_end_base','primer_type']
    dtypes= {'read_core_id': str,
    'mRNA': str,
    'init_polya_start_base': int,
    'init_polya_end_base': int,
    'primer_type': str}
    
    df = pd.read_csv(inadapter, delimiter = "\t", usecols=fields,dtype=dtypes)
    print("ok")
    
    
    df['readname'] = df['read_core_id'].str.split(",",n=1).str[0]
    #fq = SeqIO.to_dict(SeqIO.parse(gzip.open(inseq, "rt"),'fastq'))
    readnames=df['readname'].to_list
    fq={}

    with gzip.open(inseq, "rt") as handle:
        for record in tqdm(SeqIO.parse(handle, "fastq")):
            fq[record.id] = record.seq

    tqdm.pandas()

    df['read_seqs']=df.swifter.apply(lambda x: get_seq(fq, x.readname, x.primer_type in REV ), axis=1)
    print("111111___ TEST ___111111")
    print(df.head())
    print(df.columns.tolist())

    df[['polytail','additional_tail','adapter', 'dist_adapter']] = df.swifter.apply(get_three_primes_parts_row, inseq=inseq, adapt_seq=adapt_seq, axis=1 )
    df = pd.concat([df, df.apply(lambda col: get_composition(col["polytail"], "A-tail"), axis=1, result_type="expand")], axis = 1)
    df = pd.concat([df, df.apply(lambda col: get_composition(col["additional_tail"], "add-tail"), axis=1, result_type="expand")], axis = 1)
    df.drop('read_seqs', axis=1, inplace=True)
    df = df.replace(r'^\s*$', np.nan, regex=True)
    df.to_csv(out, encoding='utf-8', index=False, sep='\t', na_rep="NA")


def get_three_primes_parts_row(row, inseq, adapt_seq, debug=False):
       
    read_seq=row['read_seqs']
    primer_type=row['primer_type']

    if row['init_polya_start_base'] == row['init_polya_end_base']:
        
        return pd.Series([np.nan, np.nan, np.nan, np.nan])

    polytail=np.nan
    additional_tail=np.nan
    if primer_type in FWD:

        gene = read_seq[:row['init_polya_start_base']-1]
        polytail = read_seq[row['init_polya_start_base']-1:row['init_polya_end_base']]
        three_p_seq=read_seq[row['init_polya_end_base']:]
        #align = get_umi_alignment(three_p_seq, three_p_motive)
        if three_p_seq.startswith(adapt_seq):

            start_adapt, end_adapt=(0,len(adapt_seq)-1)
            
        else:
            result=edlib.align(adapt_seq, three_p_seq,task="path", mode='HW')

            start_adapt, end_adapt=result['locations'][0]
        additional_tail=read_seq[row['init_polya_end_base']:row['init_polya_end_base']+start_adapt]
        adapter=read_seq[row['init_polya_end_base']+start_adapt:row['init_polya_end_base']+end_adapt+1]

    elif primer_type in REV:

        gene = read_seq[:len(read_seq) -row['init_polya_end_base']]
        polytail = read_seq[len(read_seq) -row['init_polya_end_base']: len(read_seq) -row['init_polya_start_base']+1]
        three_p_seq=read_seq[len(read_seq) - row['init_polya_start_base']+1:]            
        #align = get_umi_alignment(three_p_seq, three_p_motive)
        if three_p_seq.startswith(adapt_seq):

            start_adapt, end_adapt=(0,len(adapt_seq)-1)
            
        else:
            result=edlib.align(adapt_seq, three_p_seq,task="path", mode='HW')
            start_adapt, end_adapt=result['locations'][0]
        additional_tail=read_seq[len(read_seq) - row['init_polya_start_base']+1:len(read_seq) - row['init_polya_start_base']+1+start_adapt]
        adapter=read_seq[len(read_seq) -row['init_polya_start_base']+1+start_adapt: len(read_seq) -row['init_polya_start_base']+1+end_adapt+1]
            
    else :
        raise Exception(f"Unknown primer type: {primer_type}\n{FWD}\n{REV}")
        
    dist_adapter=edlib.align(adapter, adapt_seq,task="path", mode='HW')["editDistance"]


    reconstructed_seq=gene+polytail+additional_tail+ adapter
    if debug:
        print("")
        print("###########")
        print("gene")
        print(gene)
        print("polytail")
        print(polytail)
        print("additional_tail")
        print(additional_tail)
        print("adapter")
        print(adapter)


    try:
        assert(reconstructed_seq in read_seq), "A Message"
    except AssertionError:
        print(f"#@#@#@#@#@#@#@#@ {row['readname']}, {row['primer_type']} #@#@#@#@#@#@#@#@#@#@#@#@#@")
        #print(readname, primer_type)
        print(row['mRNA'], row['init_polya_start_base'], row['init_polya_end_base'], row['primer_type'])
        print(read_seq)
        print("")

        print("gene")
        print(gene)
        print("polytail")
        print(polytail)
        print("additional_tail")
        print(additional_tail)
        print("adapter")
        print(adapter)
        return pd.Series(["assert_error", "assert_error", "assert_error", "assert_error"])



    return pd.Series([polytail, additional_tail, adapter, dist_adapter])
   
def parse_fastq(fastq):
    with open(fastq, 'r') as f:
        content = f.readlines()

        # Recreate content without lines that start with @ and +
        content = [line for line in content if not line[0] in '@+']

        # Now the lines you want are alternating, so you can make a dict
        # from key/value pairs of lists content[0::2] and content[1::2]
        data = dict(zip(content[0::2], content[1::2]))

    return data

def parse_fastqgz(fastqgz):
    with gzip.open(fastqgz, 'rt') as f:
        fastq_iterator = (l.rstrip() for l in f)
        for record in zip_longest(*[fastq_iterator] * 4):
            yield record

def get_umi_alignment(three_p_seq, subject, task="path", mode='HW'):
    wildcard=[]
    seq = subject
    for c in 'actgACTG':
        seq = seq.replace(c, "")
    wildcard = set(''.join(seq))

    result=edlib.align(subject, three_p_seq, task=task, mode=mode, additionalEqualities=equalities)
    align = edlib.getNiceAlignment(result, subject, three_p_seq)

    umi=""

    for q, t in zip(align["query_aligned"], align["target_aligned"]):
        if q not in wildcard:
            continue
        if t == "-":
            umi += "N"
        else:
            umi += t
    result["umi"] = umi

    return(result)


def get_seq(fq, readname, antisense) -> str:

    seq=fq[readname]
    if antisense:
       return str(seq.reverse_complement())
    else:
        return str(seq)

def get_composition(seq: str, col_prefix) -> pd.Series:
    """
    JP
    Get the composition and the percentage of ATGC in a DNA sequence. 
    Ignore all nucleotides other than ATGC (like N).
    """
    nucl_count = {col_prefix+'_A': 0, col_prefix+'_T': 0, col_prefix+'_G': 0, col_prefix+'_C': 0}
    nucl_perc = {col_prefix+'_pct_A': 0, col_prefix+'_pct_T': 0, col_prefix+'_pct_G': 0, col_prefix+'_pct_C': 0}


    if seq and not pd.isna(seq):
        seq = regex.sub('[^ATCG]', '', seq)
        seq_len=len(seq)
        
        if seq_len>0:
            res = Counter(seq)
            #print(res['A'])
            nucl_count = {col_prefix+'_A': res['A'], col_prefix+'_T': res['T'], col_prefix+'_G': res['G'], col_prefix+'_C': res['C']}
            nucl_perc = {col_prefix+'_pct_A': res['A']/seq_len, col_prefix+'_pct_T': res['T']/seq_len, col_prefix+'_pct_G': res['G']/seq_len, col_prefix+'_pct_C': res['C']/seq_len}

    return_val = pd.concat([pd.Series(nucl_count), pd.Series(nucl_perc)])
    return return_val



if __name__ == '__main__':
    main()