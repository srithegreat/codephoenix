__author__ = 'srikanth'

from Bio import SeqIO
def convert_tba_format():
    for i in SeqIO.parse(open("Anopheles-albimanus-STECLA_SCAFFOLDS_AalbS1.fa"),"fasta"):
        header = i.id.





if __name__ == "__main__":
    __author__ = 'Srikanth'
    parser = argparse.ArgumentParser(description='Program to parse...')
    parser.add_argument('-d', '--database', help='Input db', required=True)
    parser.add_argument('-i', '--inputfile', help='Input file', required=True)
    #parser.add_argument('-o','--output',help='Output file name', required=True)
    args = parser.parse_args()