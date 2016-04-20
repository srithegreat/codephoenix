from pyteomics import parser
import argparse
from Bio import SeqIO
import sys
def cleave_peptide(db):
    c=0
    #dict_of_sequences={}
    for i in SeqIO.parse(open(db), "fasta"):
        sequence = str(i.seq)
        peptides = parser.cleave(sequence,parser.expasy_rules['trypsin'])
        size = [j for j in peptides if (len(j)>6 and len(j)<=30)]
        #dict_of_sequences[i.id.split("|")[1]]=len(size)
	for pep in size:
		c+=1
		print str(c)+"\t"+i.id.split("|")[1]+"\t"+pep
		
	#print i.id.split("|")[1]+"\t"+i.description.split("|")[4].strip()+"\t"+str(len(i.seq))+"\t"+str(len(size))


cleave_peptide(sys.argv[1])
