#!/usr/bin/env python3
import pandas as pd
import click
import pyfastx
import sys
import edlib
import regex
from tqdm import tqdm
from colorama import Fore, Style, init
from collections import Counter
init(autoreset=True)


equalities=[("M", "A"), ("M", "C"),("R", "A"), ("R", "A"), ("W", "A"), ("W", "A"), ("S", "C"), ("S", "C"), ("Y", "C"), ("Y", "C"), 
("K", "G"), ("K", "G"), ("V", "A"), ("V", "C"), ("V", "G"), ("H", "A"), ("H", "C"), ("H", "T"), ("D", "A"), ("D", "G"), ("D", "T"),
("B", "C"), ("B", "G"), ("B", "T"), ("N", "G"), ("N", "A"), ("N", "T"), ("N", "C")]
    

FWD=["F-R", "F-UR", "F-UUR", "F-UF", "F-F", "F-N", "N-R", "N-N", "N-UR", "N-UUR", "UF-R", "UF-N", "UF-UF", "UF-UR", "UF-UUR", "UR-R", "UUR-R", "UUR-UR"]
REV=["R-F", "R-R", "R-N", "R-UF", "R-UUR", "R-UR", "N-F", "N-UF", "UF-F", "UR-F", "UR-N", "UR-UF", "UR-UR", "UR-UUR", "UUR-F", "UUR-N", "UUR-UF", "UUR-UUR"]
 
@click.command()
@click.option('-i', '--inadapter', help='Input adapter file', required=True, type=click.Path(exists=True))                         
@click.option('-s', '--inseq', help='Input read fastq file', required=True, type=click.Path(exists=True))                             
@click.option('-o', '--out', help='Output adapter information of each read', required=True)
@click.option('-v', '--verbose', is_flag=True, help="Print sequences with colors for polytail, adapter and delimiter")
@click.option('-d', '--debug', is_flag=True, help="for developping purposes, prints additional information")

def get_three_primes_parts(inadapter, inseq, out, constant_seq="CTGAC", umi_seq="NNNNNNNNNN", adapt_seq="CTGTAGGCACCATCAAT", verbose=False, debug=False):
    

    df = pd.read_csv(inadapter, delimiter = "\t")

    read_core_id_list, mRNA_list, polytail_list, additional_tail_list, umi_list, adapter_list, dist_adapter_list = [], [], [], [], [], [], []
    Count=0
    CountFwd=0
    CountRev=0
    nb_proc=0

    for readid, mRNA, polyA_start, polyA_end, primer_type in tqdm(zip(df['read_core_id'], df['mRNA'], df['init_polya_start_base'], df['init_polya_end_base'], df['primer_type']), total = df.shape[0]):
        #input("wait input\n\n")
        Count+=1

        readname = readid.split(",")[0]
        read_seq=get_seq(inseq, readname, primer_type in REV )
        print(readname, primer_type)
        print(read_seq)

        if primer_type in FWD:
            CountFwd+=1
            gene = read_seq[:polyA_start-1]
            polytail = read_seq[polyA_start-1:polyA_end]

            three_p_seq=read_seq[polyA_end:]
            #align = get_umi_alignment(three_p_seq, three_p_motive)
            result=edlib.align(adapt_seq, three_p_seq,task="path", mode='HW')
            align=edlib.getNiceAlignment(result, adapt_seq, three_p_seq)
            print(result)
            print(align)


            start_adapt, end_adapt=result['locations'][0]
            print(start_adapt, end_adapt)
            print(polyA_end)
            additional_tail=read_seq[polyA_end:polyA_end+start_adapt]
            adapter=read_seq[polyA_end+start_adapt:polyA_end+end_adapt+1]

        elif primer_type in REV:
            

            CountRev+=1
            gene = read_seq[:len(read_seq) -polyA_end]
            polytail = read_seq[len(read_seq) -polyA_end: len(read_seq) -polyA_start+1]
            three_p_seq=read_seq[len(read_seq) - polyA_start+1:]            
            #align = get_umi_alignment(three_p_seq, three_p_motive)
            result=edlib.align(adapt_seq, three_p_seq,task="path", mode='HW')
            align=edlib.getNiceAlignment(result, adapt_seq, three_p_seq)
            print(result)
            print(align)

            start_adapt, end_adapt=result['locations'][0]
            print(start_adapt, end_adapt)
            print(len(read_seq) - polyA_start+1)
            print(len(read_seq) - polyA_start+1+start_adapt)

            additional_tail=read_seq[len(read_seq) - polyA_start+1:len(read_seq) - polyA_start+1+start_adapt]
            adapter=read_seq[len(read_seq) -polyA_start+1+start_adapt: len(read_seq) -polyA_start+1+end_adapt+1]
            

       
        else :
            raise Exception(f"Unknown primer type: {primer_type}")
        
        dist_adapter=edlib.align(adapter, adapt_seq,task="path", mode='HW')["editDistance"]

        #umi = align["umi"]
        #delimiter, adapter= three_p_seq.partition(umi)[0][len(additional_tail):], three_p_seq.partition(umi)[2][:len(adapt_seq)]            
        
        #dist_delim=edlib.align(delimiter, constant_seq,task="path", mode='HW')["editDistance"] 
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


        assert(reconstructed_seq in read_seq)

        read_core_id_list.append(readid)
        mRNA_list.append(mRNA)
        polytail_list.append(polytail)
        additional_tail_list.append(additional_tail)

        adapter_list.append(adapter)
        dist_adapter_list.append(dist_adapter)        
        
        if verbose:
            print(f">{readname},-")
            print(f"{Fore.GREEN}{gene}{Fore.LIGHTYELLOW_EX}{polytail}{Fore.GREEN}{Style.RESET_ALL}{Fore.GREEN}{adapter}")
            
    #print(Counter(dist_adapter_list))
    
    res = pd.DataFrame(
        {'read_core_id':read_core_id_list,
          'gene': mRNA_list, 
          'polytail': polytail_list,
          'additional_tail' : additional_tail_list,
          'adapter': adapter_list,
          'dist_adapter_list': dist_adapter_list
          } 
    )

    print(Count, "mRNAs")
    print(CountFwd, "fwd")
    print(CountRev, "rev")
    print(f"{len(umi_list)} identified 3', with {len(set(umi_list))} unique UMIs")


    print("Adding base composition")
    res = pd.concat([res, res.apply(lambda col: get_composition(col["polytail"], "A-tail"), axis=1, result_type="expand")], axis = 1)
    res = pd.concat([res, res.apply(lambda col: get_composition(col["additional_tail"], "add-tail"), axis=1, result_type="expand")], axis = 1)


    res.to_csv(out, encoding='utf-8', index=False)

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


def get_seq(fastq, readname, antisense=False) -> str:
    
    fq=pyfastx.Fastq(fastq)

    seq = fq[readname].seq if not antisense else fq[readname].antisense
    return seq

def get_composition(seq: str, col_prefix) -> pd.Series:
    """
    JP
    Get the composition and the percentage of ATGC in a DNA sequence. 
    Ignore all nucleotides other than ATGC (like N).
    """
    nucl_count = {col_prefix+'_A': 0, col_prefix+'_T': 0, col_prefix+'_G': 0, col_prefix+'_C': 0}
    nucl_perc = {col_prefix+'_pct_A': 0, col_prefix+'_pct_T': 0, col_prefix+'_pct_G': 0, col_prefix+'_pct_C': 0}

    seq = regex.sub('[^ATCG]', '', seq)
    seq_len=len(seq)
    if seq and not pd.isna(seq):
        
        if seq_len>0:
            res = Counter(seq)
            #print(res['A'])
            nucl_count = {col_prefix+'_A': res['A'], col_prefix+'_T': res['T'], col_prefix+'_G': res['G'], col_prefix+'_C': res['C']}
            nucl_perc = {col_prefix+'_pct_A': res['A']/seq_len, col_prefix+'_pct_T': res['T']/seq_len, col_prefix+'_pct_G': res['G']/seq_len, col_prefix+'_pct_C': res['C']/seq_len}

    return_val = pd.concat([pd.Series(nucl_count), pd.Series(nucl_perc)])
    return return_val



if __name__ == '__main__':
    get_three_primes_parts()


