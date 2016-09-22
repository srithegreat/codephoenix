__author__ = 'Srikanth'
from pyteomics import parser
import argparse
from Bio import SeqIO
import math
def cleave_peptide(db,inp):
    dict_of_sequences={}
    for i in SeqIO.parse(open(db), "fasta"):
        sequence = str(i.seq)
        peptides = parser.cleave(sequence,parser.expasy_rules['trypsin'])
        size = [j for j in peptides if (len(j)>6 and len(j)<=35)]
        dict_of_sequences[i.id.split("|")[1]]=len(size)


    summed = {}
    for j1 in open(inp):
        if not j1.startswith("Sequence"):
            peptide, accession, intensity = j1.strip().split("\t")

            for k in accession.split(";"):
                if k in summed:
                    summed[k] += float(intensity)
                else:
                    summed[k]=0
                    summed[k] += float(intensity)

    for each in summed:
        lengths = dict_of_sequences.get(each,0)
        if lengths==0:
            ibaq=0
        else:
            ibaq = float(summed[each])/lengths
        try:
            print "\t".join([str(each), str(ibaq),str(math.log10(ibaq))])
        except:
            print "\t".join([str(each), str(ibaq),"0"])

#cleave_peptide("/home/iob/Dropbox/Sharing/6_Tissue_Transcriptome/Sequences/Gencode_Hs_V22_21816.fasta","/home/iob/Dropbox/Sharing/6_Tissue_Transcriptome/Transcriptomic/Data.txt")#requires are fasta genome file and data file with three columns

if __name__ == "__main__":
    __author__ = 'Srikanth'
    parser = argparse.ArgumentParser(description='Program to get iBAQ values, print to Std out...')
    parser.add_argument('-d', '--database', help='Input db', required=True)
    parser.add_argument('-i', '--inputfile', help='Input file', required=True)
    #parser.add_argument('-o','--output',help='Output file name', required=True)
    args = parser.parse_args()
    cleave_peptide(args.database,args.inputfile)
