__author__ = 'iob'
from Bio import SeqIO
import re, argparse


def getseq(genomefasta):
    """ Returns the dictionary of genome fasta"""
    genomedict = {}
    for i in SeqIO.parse(open(genomefasta), "fasta"):
        genomedict[i.id] = str(i.seq)
    return genomedict


def read_gff(gff):
    """Program to read a gff and create dictionary of exons from a transcript"""
    genome = getseq(args.genome)
    dictoftranscripts = {}
    for k in open(gff):
        if not k.startswith("#"):
            lines = k.strip().split("\t")
            if lines[2] == "exon":
                strand = lines[6]
                chromosome = lines[0]
                start = lines[3]
                end = lines[4]
                transcriptid = re.search("Parent=transcript:(.*)", lines[8]).group(1)
                if transcriptid + "#" + chromosome in dictoftranscripts:
                    dictoftranscripts[transcriptid + "#" + chromosome].extend([start, end])
                else:
                    dictoftranscripts[transcriptid + "#" + chromosome] = []
                    dictoftranscripts[transcriptid + "#" + chromosome].extend([start, end])

    for key, value in dictoftranscripts.iteritems():
        value.sort()
        print value
        for coord1 in value:

            for coord2 in value[1:]:
                #print coord1, coord2
                if int(coord1) != int(value[-1]) and value.index(coord2) != value.index(coord1)+1 and value.index(coord2) > value.index(coord1):

                    exon1_start = int(coord1)
                    exon1_end = int(coord2)
                    #print exon1_start, exon1_end
                    #print key.split("#")[1]
                    #print value.index(coord1), value.index(coord2)
                    exon_seq = genome.get(key.split("#")[1],"NA")

                    if exon_seq != "NA":
                        sequence_exon = exon_seq[exon1_start:exon1_end+1]
                        #print exon1_start, exon1_end, sequence_exon
                        for start, end, strand, frame, pro in translate(sequence_exon):
                            junction =
                            print start, end, strand, frame, pro


def translate(sequence):
    min_protein_length = 6
    from Bio.Seq import Seq
    s = Seq(sequence)
    seq_len = len(s)
    answer = []
    for strand, nucleotide in [(+1, s)]:
            for frame in range(3):
                if len(nucleotide[frame:])% 3 != 0:
                    trans=str((nucleotide[frame:]+"N"*(3-len(nucleotide[frame:])%3)).translate())
                else:
                    trans = str(nucleotide[frame:].translate())
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





if __name__ == "__main__":
    __author__ = 'Srikanth'
    parser = argparse.ArgumentParser(description='Program to create exon shuffle database...')
    parser.add_argument('-g', '--genome', help='genome fasta', required=True)
    parser.add_argument('-f', '--gff', help='gff file', required=True)
    # parser.add_argument('-o','--output',help='Output file name', required=True)
    args = parser.parse_args()
    read_gff(args.gff)
