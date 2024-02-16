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

# FLEPSEQ2

equalities=[("M", "A"), ("M", "C"),("R", "A"), ("R", "A"), ("W", "A"), ("W", "A"), ("S", "C"), ("S", "C"), ("Y", "C"), ("Y", "C"), 
("K", "G"), ("K", "G"), ("V", "A"), ("V", "C"), ("V", "G"), ("H", "A"), ("H", "C"), ("H", "T"), ("D", "A"), ("D", "G"), ("D", "T"),
("B", "C"), ("B", "G"), ("B", "T"), ("N", "G"), ("N", "A"), ("N", "T"), ("N", "C")]
    

FWD=("F-R", "F-UR", "F-UUR", "F-UF", "F-F", "F-N", "N-R", "N-N", "N-UR", "N-UUR", "UF-R", "UF-N", "UF-UF", "UF-UR", "UF-UUR", "UR-R", "UUR-R", "UUR-UR")
REV=("R-F", "R-R", "R-N", "R-UF", "R-UUR", "R-UR", "N-F", "N-UF", "UF-F", "UR-F", "UR-N", "UR-UF", "UR-UR", "UR-UUR", "UUR-F", "UUR-N", "UUR-UF", "UUR-UUR")
 
@click.command()
@click.option('--inadapter', help='Input adapter file', required=True, type=click.Path(exists=True))                         
@click.option('--inseq', help='Input read fastq file', required=True, type=click.Path(exists=True))                             
@click.option('--out', help='Output adapter information of each read', required=True)
@click.option('--max_tail', help="Max 3' sequence length (polyA + addTail) to consider (avoids too intensive alignment calculations) ", required=True, type=int)
@click.option('--debug', is_flag=True, help="for developping purposes, prints additional information")

def main(inadapter, inseq, out, max_tail, debug, constant_seq="CTGAC", umi_seq="NNNNNNNNNN", adapt_seq="CTGTAGGCACCATCAAT"):
    fields = ['read_core_id', 'chr', 'read_exon_total_num','mRNA', "mRNA_start", "mRNA_end", 'mRNA_intron_num', 'retention_introns',
              'polya_start_raw', 'polya_end_raw','polya_start_base', 'polya_end_base', 
              'init_polya_start_base', 'init_polya_end_base','primer_type', 'polya_length', 'init_polya_length', 'type']
    dtypes= {'read_core_id': str,
             'chr': str,
             'read_exon_total_num': int,
             'mRNA': str,
             'mRNA_start': int,
             'mRNA_end' : int,
             'mRNA_intron_num': int,
             'retention_introns': str,
             'polya_start_base': int,
             'polya_end_base':int,
             'init_polya_start_base': int,
             'init_polya_end_base': int,
             'primer_type': str,
             'polya_length': float,
             'init_polya_length': float,
             'type': str}
    
        
    log_dict={}
    df = pd.read_csv(inadapter, delimiter = "\t", usecols=fields,dtype=dtypes)
    log_dict["nb_reads_input"]=len(df.index)     
    
    df['readname'] = df['read_core_id'].str.split(",",n=1).str[0]
    #fq = SeqIO.to_dict(SeqIO.parse(gzip.open(inseq, "rt"),'fastq'))
    readnames=df['readname'].to_list
    fq={}

    with gzip.open(inseq, "rt") as handle:
        for record in tqdm(SeqIO.parse(handle, "fastq")):
            fq[record.id] = record.seq

    tqdm.pandas()
    df['read_seqs']=df.swifter.apply(lambda x: get_seq(fq, x.readname, x.primer_type in REV ), axis=1)
    df['sense']= np.where(df["primer_type"].isin(FWD), "FWD", "REV")
    # This step is to avoid inconsistencies when calculating boundaries (polya, additional tail...)
    df['init_polya_end_base'] = np.where((df['init_polya_start_base'] == df['init_polya_end_base']) , df['init_polya_end_base']-1, df['init_polya_end_base'])
    df[['polytail','additional_tail','adapter', 'dist_adapter','coords_in_read', 'comment']] = df.apply(get_three_primes_parts_row, max_tail=max_tail, adapt_seq=adapt_seq, axis=1, debug=debug)
    #df[['polytail','additional_tail','adapter', 'dist_adapter','coords_in_read', 'comment']] = df.swifter.apply(get_three_primes_parts_row, adapt_seq=adapt_seq, axis=1, debug=debug)
    df = pd.concat([df, df.apply(lambda col: get_composition(col["polytail"], "A_tail"), axis=1, result_type="expand")], axis = 1)
    df = pd.concat([df, df.apply(lambda col: get_composition(col["additional_tail"], "add_tail"), axis=1, result_type="expand")], axis = 1)
    df.drop('read_seqs', axis=1, inplace=True)
    df = df.replace(r'^\s*$', np.nan, regex=True)
    df.to_csv(out, encoding='utf-8', index=False, sep='\t', na_rep="NA")
    
    #Making log
    log_dict["nb_reads_output"]=len(df.index)
    log_dict["nb_unid_polyA"]=df['polytail'].isna().sum()
    log_dict["nb_unid_add"]=df['additional_tail'].isna().sum()
    log_dict["nb_unid_adapter"]=df['adapter'].isna().sum()
    log_dict["nb_reads_fwd"]=df.sense.value_counts()['FWD']
    log_dict["nb_reads_rev"]=df.sense.value_counts()['REV']

    with open(out+'.log', 'w') as outlog:
        print(log_dict, file=outlog)

def get_three_primes_parts_row(row, max_tail, adapt_seq, debug):
    read_seq=row['read_seqs']
    primer_type=row['primer_type']
    comment="everything OK"

    polytail=np.nan
    additional_tail=np.nan
    try:
        if primer_type in FWD:            
            gene_start_in_read=0
            gene_end_in_read=row['init_polya_start_base']-1
            gene = read_seq[:gene_end_in_read]
            polytail_start_in_read=row['init_polya_start_base']-1
            polytail_end_in_read=row['init_polya_end_base']
            polytail = read_seq[polytail_start_in_read:polytail_end_in_read]
            three_p_seq_start_in_read=row['init_polya_end_base']
            three_p_seq_end_in_read=len(read_seq)
            three_p_seq=read_seq[three_p_seq_start_in_read:]
            if len(three_p_seq)>max_tail:
                three_p_seq_end_in_read=max_tail
                three_p_seq=three_p_seq[:three_p_seq_end_in_read]
                comment="3' seq detected >max_tail, shortened"   
                
            #align = get_umi_alignment(three_p_seq, three_p_motive)
            if three_p_seq.startswith(adapt_seq):

                start_adapt, end_adapt=(0,len(adapt_seq)-1)
                
            else:
                result=edlib.align(adapt_seq, three_p_seq,task="path", mode='HW')

                start_adapt, end_adapt=result['locations'][0]
            add_tail_start_in_read=row['init_polya_end_base']
            add_tail_end_in_read=row['init_polya_end_base']+start_adapt
            additional_tail=read_seq[add_tail_start_in_read:add_tail_end_in_read]
            adapt_start_in_read = row['init_polya_end_base']+start_adapt
            adapt_end_in_read = row['init_polya_end_base']+end_adapt+1
            adapter=read_seq[adapt_start_in_read:adapt_end_in_read]

        elif primer_type in REV:
            gene_start_in_read=0
            gene_end_in_read=len(read_seq) -row['init_polya_end_base']
            gene = read_seq[:gene_end_in_read]
            polytail_start_in_read=len(read_seq) -row['init_polya_end_base']
            polytail_end_in_read=len(read_seq) -row['init_polya_start_base']+1
            polytail = read_seq[polytail_start_in_read: polytail_end_in_read]
            three_p_seq_start_in_read=len(read_seq) - row['init_polya_start_base']+1
            three_p_seq_end_in_read=len(read_seq)
            three_p_seq=read_seq[three_p_seq_start_in_read:]   
            if len(three_p_seq)>max_tail:
                three_p_seq_end_in_read=max_tail
                three_p_seq=three_p_seq[:three_p_seq_end_in_read]
                comment="3' seq detected >max_tail, shortened"         
            #align = get_umi_alignment(three_p_seq, three_p_motive)
            if three_p_seq.startswith(adapt_seq):

                start_adapt, end_adapt=(0,len(adapt_seq)-1)
                
            else:
                result=edlib.align(adapt_seq, three_p_seq,task="path", mode='HW')
                start_adapt, end_adapt=result['locations'][0]
            #print("adapt start end")
            #print(start_adapt, end_adapt)

            add_tail_start_in_read=len(read_seq) - row['init_polya_start_base']+1
            add_tail_end_in_read =len(read_seq) - row['init_polya_start_base']+1+start_adapt
            additional_tail=read_seq[add_tail_start_in_read:add_tail_end_in_read]
            adapt_start_in_read=len(read_seq) -row['init_polya_start_base']+1+start_adapt
            adapt_end_in_read=len(read_seq) -row['init_polya_start_base']+1+end_adapt+1
            adapter=read_seq[adapt_start_in_read:adapt_end_in_read]
                
        else :
            raise Exception(f"Unknown primer type: {primer_type}\n{FWD}\n{REV}")
    except Exception as e:
        return pd.Series([np.nan, np.nan, np.nan, np.nan,"no coords", f"{e} occurred"])
        
    dist_adapter=edlib.align(adapter, adapt_seq,task="path", mode='HW')["editDistance"]
    reconstructed_seq=gene+polytail+additional_tail+ adapter

    reconstructed_coords_in_read=f"{gene_start_in_read}:{gene_end_in_read}:{polytail_start_in_read}:{polytail_end_in_read}:{add_tail_start_in_read}:{add_tail_end_in_read}:{adapt_start_in_read}:{adapt_end_in_read}"
    
    # print(debug)
    if debug:
        print("")
        print("###########")
        print(row)
        print("------------------")
        print(read_seq)
        print(primer_type)
        print("gene")
        print(gene)
        print("polytail")
        print(polytail)
        print("additional_tail")
        print(additional_tail)
        print("adapter")
        print(adapter)
        input("press a key")


    try:
        assert(reconstructed_seq in read_seq), "A Message"
        
    except AssertionError:
        sys.exit("Assertion error for reconstructed sequence")

    return pd.Series([polytail, additional_tail, adapter, dist_adapter,reconstructed_coords_in_read, comment])
   
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
    nucl_count = {col_prefix+'_A': np.nan, col_prefix+'_T': np.nan, col_prefix+'_G': np.nan, col_prefix+'_C': np.nan}
    nucl_perc = {col_prefix+'_pct_A': np.nan, col_prefix+'_pct_T': np.nan, col_prefix+'_pct_G': np.nan, col_prefix+'_pct_C': np.nan}

    if seq and not pd.isna(seq):

        seq = regex.sub('[^ATCG]', '', seq)
        seq_len=len(seq)  
        if seq_len>0:
            res = Counter(seq)
            #print(res['A'])
            nucl_count = {col_prefix+'_A': res['A'], col_prefix+'_T': res['T'], col_prefix+'_G': res['G'], col_prefix+'_C': res['C']}
            nucl_perc = {col_prefix+'_pct_A': res['A']*100/seq_len, col_prefix+'_pct_T': res['T']*100/seq_len, col_prefix+'_pct_G': res['G']*100/seq_len, col_prefix+'_pct_C': res['C']*100/seq_len}
    
    return_val = pd.concat([pd.Series(nucl_count), pd.Series(nucl_perc)])
    
    return return_val


if __name__ == '__main__':
    main()