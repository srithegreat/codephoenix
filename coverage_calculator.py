# -*- coding: utf-8 -*-
"""
Created on Wed Nov 20 13:40:09 2013

@author: srikanth
Calculates Coverage for any given peptides
"""
from __future__ import division
from Bio import SeqIO
import argparse
def coverage_calculator(inp):
    nomatch=open("nomatching.txt",'w')
    dictofseq={}#contains accession and corresponding sequences
    for k in open(inp):
        sequence, pep = k.strip().split("\t")
        positions=[]
        sequencelength = len(sequence)
        for every in pep.split(";"):
            if sequence.find(every)!=-1:
                positions.extend(range(sequence.index(every),sequence.index(every)+len(every)))
                totalaminocovered=len(set(positions))
            else:
                nomatch.write(k)
        #create unique set of positions
            #print sequencelength,totalaminocovered
        print "%s\t%.2f"%(sequence,(totalaminocovered/sequencelength)*100) #calculate coverage to two decimal places

    
    
if __name__ == "__main__":
    __author__ = 'Srikanth'
    parser = argparse.ArgumentParser(description='Program to calculate coverage')
    #parser.add_argument('-d', '--database', help='Input db (fasta)', required=True)
    parser.add_argument('-i', '--inputfile', help='Input file with peptides', required=True)
    #parser.add_argument('-o','--output',help='Output file name', required=True)
    args = parser.parse_args()
    coverage_calculator(args.inputfile)

       
        
