# -*- coding: utf-8 -*-
"""


@author: srikanth
Requires Biopython, TopHat, Cufflinks and CPAT in path
Replace the CPAT hexamer and logit files to actual paths before running
Works with merged.gtf from cuffmerge
"""

import pylab as P
import re,os,random,sys
from Bio import SeqIO


def main(genome, gtf, output, fpkm):
	
	
    print "Working with ", gtf
    directory = "temp"

    if not os.path.exists(directory):
        os.makedirs(directory)

    filter_gtf(gtf)
    print "Filtered GTF for class codes, creating fasta"
    temp = "temp/temp_" + str(random.randint(1,100)) + ".fasta"
    #os.system("gffread %s -U -g %s -w %s"%("filtered_gtffile.gtf",genome,temp))#single exon -U
    os.system("gffread %s -U -g %s -w %s"%("temp/filtered_gtffile.gtf", genome, temp))#Exclude single exon
    print "Done fasta, Filtering single exons and length>200 now...."
    filter_fasta(temp, output, fpkm)
    print 'Done........'
    

def filter_gtf(gtf):
    qualified={}
    print "filtering class codes of interest"
    filtered=open("temp/filtered_gtffile.gtf", "w")
    for lines in open(gtf):
        transcriptid= re.search("transcript_id\s+\"(.*?)\"",lines.strip().split("\t")[8]).group(1)
        geneid= re.search("gene_id\s+\"(.*?)\"",lines.strip().split("\t")[8]).group(1)
        classcode= re.search("class_code\s+\"(.*?)\"",lines.strip().split("\t")[8]).group(1)
        if classcode=="u": #or classcode=="x" or classcode == "i":# integenic=u,antisense=x,intronic=i
            filtered.write(lines)

    filtered.close()

    return 0
    
def filter_fasta(fname, output, fpkm):
    filtered_fasta=open("temp/filtered_transcripts.fasta","w")
    length=200
    records=0
    tc_seq={}
    for each in SeqIO.parse(open(fname),"fasta"):
        if len(each.seq)>length:
            records+=1
            filtered_fasta.write(">"+each.id+"\n"+str(each.seq)+"\n")
            tc_seq[each.id]= str(each.seq)
    filtered_fasta.close()
    print "Written records", records
    cpat_calc(tc_seq, output, fpkm)
    
    
def cpat_calc(mapped, output, fpkm):
    allids={}
    header= open(fpkm).readline().strip()
    for e in open(fpkm):
        allids[e.strip().split("\t")[0]]= e.strip().split("\t")[1:]
    print "calculating CPAT and filtering..."
    total= 0
    #out1=open("final_list_sorted"+userinp[:-4]+".txt","w")
    os.system("cpat.py -g temp/filtered_transcripts.fasta -o %s -x %s/Human_Hexamer.tsv -d %s/Human_logitModel.RData "%("temp/tempout","/home/iob/Desktop/CPAT-1.2.2","/home/iob/Desktop/CPAT-1.2.2"))#Replace pathtoCPAThexamerfile,pathtoCPATlogitmodel with actual paths before running
    prob_transcripts_multiple=open("Putative_lncRNAs_Multiple_"+output+".fasta",'w')
    filtered_transcripts=open("Putative_TUCPs.fasta","w")
    fpkm_values=open("FPKMs_putative_lncRNAs_"+output+".txt", "w")
    fpkm_values.write(header+"\n")
    filterid=[]
    for each in open("temp/tempout"):
        try:
            if float(each.strip().split("\t")[5]) < 0.364:#decided based on datasets

                filterid.append(each.strip().split("\t")[0])
        except:
            pass

    for every in mapped:
        if every in filterid:

            allfpkm=max([float(g) for g in allids.get(every, "0")])
            if allfpkm>=0.1:
                total+= 1
                fpkm_values.write(every+"\t"+"\t".join(allids[every])+"\n")
                prob_transcripts_multiple.write(">"+every+"\n"+mapped[every]+"\n")
        else:
            filtered_transcripts.write(">"+every+"\n"+mapped.get(every,"0")+"\n")



    print "Total putative lncRNAs ",total

    
def getintergenic_only():
    os.system("intersectBed -a temp/filtered_gtffile.gtf -b gencode.v22.annotation.gtf -v>intergenic_only_novel_lncs.gtf")
    for k in open("intergenic_only_novel_lncs.gtf"):pass




if __name__ == "__main__":

    import argparse, os
    os.system("mkdir temp")

    __author__ = 'Srikanth'

    parser = argparse.ArgumentParser(description='Program to Identify novel lncRNAs')
    parser.add_argument('-g','--genomefasta', help='Genome fasta sequence',required=True)
    parser.add_argument('-r','--gtf', help='Merged GTF file with class codes',required=True)
    parser.add_argument('-o','--output',help='Output file name extension', required=True)
    parser.add_argument("-f", '--fpkm', help="FPKM table matrix of isoforms (isoforms.fpkm_table)", required=True)
    #parser.add_argument('-f','--fpkm',help='Output fpkm values as tabular (default False)', choices=["True", "False"], type=bool)
    args = parser.parse_args()
    main(args.genomefasta, args.gtf, args.output, args.fpkm)
