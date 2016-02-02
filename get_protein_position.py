

from Bio import SeqIO
import re, argparse
def read_fasta_get_pos(args):
    dict_of_seq={}
    for i in SeqIO.parse(open(args.database),"fasta"):

        dict_of_seq[i.id.split("|")[1]] = str(i.seq)
        #print i.id.split("|")[1]

    c=0
    out=open(args.output,"w")
    for line in open(args.inputfile):
        pep, acc = line.strip().split("\t")
        sequences = [dict_of_seq.get(s).index(pep.strip()) for s in acc.split(";")]



        out.write(pep.strip() + "\t" + ";".join(sequences) + "\n")
        c+=1
        if c%50==0:
            print "done", c




if __name__=="__main__":

    __author__ = 'Srikanth'
    parser = argparse.ArgumentParser(description='Program to get peptide position in protein.')
    parser.add_argument('-d','--database', help='Input RefSeq database name',required=True)
    parser.add_argument('-i','--inputfile', help='Input file (pep,acc)',required=True)
    parser.add_argument('-o','--output',help='Output file name', required=True)
    args = parser.parse_args()
    read_fasta_get_pos(args)