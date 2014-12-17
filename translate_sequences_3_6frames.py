# -*- coding: utf-8 -*-
"""
Created on Thu Aug 14 15:52:31 2014

@author: srikanth
Program to translate in three or six frames
Requires BioPython
"""
from Bio import SeqIO 
import re
import sys
table = 1
min_pro_len = 7
length_orf=20

def orfs_by_trans(seq, trans_table, min_protein_length,transframes):
    answer = []
    seq_len = len(seq)
    if transframes==6:
        for strand, nucleotide in [(+1, seq), (-1, seq.reverse_complement())]:
            for frame in range(3):
                if len(nucleotide[frame:])%3!=0:
                    trans=str((nucleotide[frame:]+"N"*(3-len(nucleotide[frame:])%3)).translate(trans_table))
                else:
                    trans = str(nucleotide[frame:].translate(trans_table))
                trans_len = len(trans)
                aa_start = 0
                aa_end = 0
                while aa_start < trans_len:
                    aa_end = trans.find("*", aa_start)
                    if aa_end == -1:
                        aa_end = trans_len
                    if aa_end-aa_start >= min_protein_length:
                        if strand == 1:
                            start = frame+aa_start*3
                            end = min(seq_len,frame+aa_end*3+3)
                        else:
                            start = seq_len-frame-aa_end*3-3
                            end = seq_len-frame-aa_start*3                        
                        answer.append((start, end, strand,frame,trans[aa_start:aa_end]))
                    aa_start = aa_end+1
    elif transframes==3:
        for strand, nucleotide in [(+1, seq)]:
            for frame in range(3):
                if len(nucleotide[frame:])%3!=0:
                    trans=str((nucleotide[frame:]+"N"*(3-len(nucleotide[frame:])%3)).translate(trans_table))
                else:
                    trans = str(nucleotide[frame:].translate(trans_table))
                trans_len = len(trans)
                aa_start = 0
                aa_end = 0
                while aa_start < trans_len:
                    aa_end = trans.find("*", aa_start)
                    if aa_end == -1:
                        aa_end = trans_len
                    if aa_end-aa_start >= min_protein_length:
                        if strand == 1:
                            start = frame+aa_start*3
                            end = min(seq_len,frame+aa_end*3+3)
                        else:
                            start = seq_len-frame-aa_end*3-3
                            end = seq_len-frame-aa_start*3                        
                        answer.append((start, end, strand,frame,trans[aa_start:aa_end]))
                    aa_start = aa_end+1
    
    return answer
#print len(sys.argv)
if len(sys.argv)!=5:
    print "Usage: python translate_sequences.py input.fasta output.fasta prefixheader #offrames"
    exit()
else:
    cnt=0
    output=open(sys.argv[2],"w")
    indexfile=open(sys.argv[2]+".index","w")
    prefix=sys.argv[3]
    frames=sys.argv[4]
    for record in SeqIO.parse(open(sys.argv[1]),"fasta"):
        orf_list = orfs_by_trans(record.seq, table, min_pro_len,int(frames))
        header=record.description
        for start, end, strand, frame, pro in orf_list:
           #print pro
           if "X" not in pro:
               #orfs=re.findall("M.*",pro)
               #for each in orfs:
               if len(pro)>=min_pro_len:
                   indexfile.write(prefix+str(cnt)+"\t"+header+"|"+str(strand)+"|"+str(start)+"-"+str(end)+":"+str(frame+1)+"\n")
                   output.write( ">"+prefix+str(cnt)+"\n"+pro+"\n")
                   cnt+=1
print "Translated to %i sequences"%cnt