#!/usr/bin/env python
# - Helene Zuber <helene.zuber@ibmp-cnrs.unistra.fr>

#################################################################################
# GOAL: Extract the portion of sequence that doesn't map the reference sequences
#1) Parse bam files extract read attributes
#2) Extract from the sequence the unmapped portion based on read.query_alignment_start or read.query_alignment_end
#3)  Clip out the adaptor and reverse complement sequence for reverse sequence
#4) Analyse composition of the tail
#################################################################################

import re
import regex
import gzip
import os, sys
import edlib
import pysam
from Bio.Seq import Seq
from collections import Counter


###################################################################################
# Main function
def main() :
    
    Bamfile_name = sys.argv[1] #read 1 fastq files
    Outfile_name =  sys.argv[2] #directory where deduplicated read 1 files will be saved

    #create an index file for bam 
    pysam.index(Bamfile_name)
    Bamfile = pysam.AlignmentFile(Bamfile_name, "rb")
    Outfile = gzip.open(Outfile_name,"wt")
        #Print input and output files    
    print(Bamfile_name)
    print(Outfile_name)
    print ("Total\tForward\tReverse\tA_onlyA\tB_onlyU\tC_onlyG\tD_onlyC\tE_polyAandU\tF_polyAandG\tG_polyAandC\tH_otherU\tI_otherG\tJ_otherC\tK_unclassified")
    
    Outfile.write("Read_id\tFlag\tReference_name\tReference_start\tReference_end\tSense\tTail\tlen_tail\tAddedNt\tlen_added\taddedNt_pctAt\taddedNt_pctC\taddedNt_pctG\taddedNt_pctU\tTag\tRead_id_custom\n")
    # Define string to be searched
    
    
    ## Search added nucleotide in the tail
    a_SearchString_onlyA = '(^A+$)'
    b_SearchString_onlyU = '(^T+$)'
    c_SearchString_onlyG = '(^G+$)'
    d_SearchString_onlyC = '(^C+$)'
    
    e_SearchStringpolyAandTail = '(^A{3,})(.*$)'

    
    ## Define counters
    Count=0
    CountFor=0
    CountRev=0
    
    
    Counta = 0
    Countb = 0
    Countc = 0
    Countd = 0
    
    Counte = 0

    Countk = 0

    adapt_seq=str(Seq("CTGTAGGCACCATCAAT"))
    adapt_seq_rc=str(Seq("CTGTAGGCACCATCAAT").reverse_complement())


    for read in Bamfile:
        #print("##################################")
        Pass = 'FALSE'
        Tail=None
        AddedNt =None
        Tag =None
        Count+=1
        Read_id = read.query_name
        Query_alignment_start = read.query_alignment_start
        Query_alignment_end = read.query_alignment_end
        Query_length = read.query_length
        Flag = read.flag
        if Flag is None:
            print(read)
        Reference_name = read.reference_name
        Reference_start = read.reference_start
        Reference_end = read.reference_end
        Query_sequence = read.query_sequence
        
        Left_Unmapped_seq = Query_sequence[:Query_alignment_start]
        Right_Unmapped_seq = Query_sequence[Query_alignment_end:]
        # print("LEFT UNMAPPED")
        # print(Left_Unmapped_seq)
        # print("QUERY SEQ")
        # print(Query_sequence)
        # print("RIGHT UNMAPPED")
        # print(Right_Unmapped_seq)


        # print("ADAPT Sequences to look for")
        # print(adapt_seq)
        # print(adapt_seq_rc)
        

        result=edlib.align(adapt_seq, Right_Unmapped_seq,task="path", mode='HW')
        #print(result)
        dist_adapter=result["editDistance"]
        #nice = edlib.getNiceAlignment(result, adapt_seq, Right_Unmapped_seq)
        result_rc=edlib.align(adapt_seq_rc, Left_Unmapped_seq,task="path", mode='HW')
        #print(result_rc)
        #print("ok")
        dist_adapter_rc=result_rc["editDistance"]
        #nice_rc = edlib.getNiceAlignment(result_rc, adapt_seq_rc, Left_Unmapped_seq)
        nb_fwd=0
        nb_rev=0
        if dist_adapter<dist_adapter_rc:
            #print("fwd!!!")
            Sense="FWD"
            nb_fwd+=1
            Tailseq = Right_Unmapped_seq[:result["locations"][0][0]]
            #print(Tail)
        else:
            #print("rev")
            Sense="REV"
            nb_rev+=1
            Tailseq=str(Seq(Left_Unmapped_seq[result_rc["locations"][0][1]+1:]).reverse_complement())
            #print(Tail)

        
        #print("-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-")
        
        ## Function
        if regex.findall(a_SearchString_onlyA, Tailseq):
            Result=regex.search(a_SearchString_onlyA, Tailseq)
            Tail=Result.group(1) 
            Tag = 'A_onlyA'
            Counta += 1
        elif regex.findall(b_SearchString_onlyU, Tailseq):
            Result=regex.search(b_SearchString_onlyU, Tailseq)
            AddedNt=Result.group(1) 
            Tag = 'B_onlyU'
            Countb += 1
        elif regex.findall(c_SearchString_onlyG, Tailseq):
            Result=regex.search(c_SearchString_onlyG, Tailseq)
            AddedNt=Result.group(1) 
            Tag = 'C_onlyG'
            Countc += 1
        elif regex.findall(d_SearchString_onlyC, Tailseq):
            Result=regex.search(d_SearchString_onlyC, Tailseq)
            AddedNt=Result.group(1) 
            Tag = 'D_onlyC'
            Countd += 1
            
        elif regex.findall(e_SearchStringpolyAandTail, Tailseq):
            Result=regex.search(e_SearchStringpolyAandTail, Tailseq)
            Tail=Result.group(1)
            AddedNt=Result.group(2) 
            Tag = 'E_polyAandTail'
            Counte += 1

        else :
            AddedNt=None if Tailseq == "" else Tailseq
            Tag = 'k_unclassified'
            Countk += 1
        
        if AddedNt is None:
            len_added=0
            addedNt_pctA, addedNt_pctC, addedNt_pctG, addedNt_pctU=0,0,0,0
        else:
            len_added=len(AddedNt)
            addedNt_pctA=Counter(AddedNt)['A']*100/len_added
            addedNt_pctC=Counter(AddedNt)['C']*100/len_added
            addedNt_pctG=Counter(AddedNt)['G']*100/len_added
            addedNt_pctU=Counter(AddedNt)['T']*100/len_added
        if Tail is None:
            len_tail=0
        else:
            len_tail=len(Tail)
        Read_id_custom= f"{Read_id},{Reference_name},{Reference_start}"
        
        #print(f"{Read_id}\t{Flag}\t{Reference_name}\t{Reference_start}\t{Reference_end}\t{Sense}\t{Tail}\t{len_tail}\t{AddedNt}\t{len_added}\t{addedNt_pctA}\t{addedNt_pctC}\t{addedNt_pctG}\t{addedNt_pctU}\t{Tag}\t{Read_id_custom}\n")
        Outfile.write(f"{Read_id}\t{Flag}\t{Reference_name}\t{Reference_start}\t{Reference_end}\t{Sense}\t{Tail}\t{len_tail}\t{AddedNt}\t{len_added}\t{addedNt_pctA}\t{addedNt_pctC}\t{addedNt_pctG}\t{addedNt_pctU}\t{Tag}\t{Read_id_custom}\n")
            
    
    Bamfile.close()
    Outfile.close()
    #print ("%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i" % (Count, CountFor, CountRev, CountA, CountB, CountC, CountD, CountE, CountF, CountG, CountH, CountI, CountJ, CountK))

    
if __name__ == "__main__" :
    main()


