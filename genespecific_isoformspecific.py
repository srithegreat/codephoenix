# -*- coding: utf-8 -*-
"""
Created on Thu Jan 15 15:20:53 2015

@author: srikanth
"""
from Bio import SeqIO
import re
import sys
refseq="/media/srikanth/My Passport/backup_system/Downloads/refseq65/Human_RefSeq65_withGS_061214.fasta"
def genespecific_isoformspecific():
    #multipleisogenes=count_isoforms_pergene(refseq)
    print "Peptide\tGene Symbol\tLength\tGenespecific\tIsoformSpecific\tShared\taccessions"
    for each in open("/home/srikanth/Downloads/LRKK2_dawson/temp_matched"):
        #print each
        peptide,accession=each.strip().split("\t")
        
        gene=[c.strip().split("|")[4].split("$")[0] for c in accession.split(";") if not c.strip().split("|")[1].startswith("c")]
        peplength=len(peptide)
        if len(set(gene))==1:
            genespec=True
            
            #if len(accession.split(";"))==1 and (gene[0] in multipleisogenes):
             #   isoformspec=True
              #  shared=False
            #else:
             #   isoformspec=False
              #  shared=False
        #elif len(set(gene))>1:
         #   genespec=False
          #  shared=True
            
        #nacc=[]
        #for i in accession.split(";"):
         #   if i.split("|")[1].startswith("c"):
          #      nacc.append(i.split("|")[1])
           # else:
            #    nacc.append(i.strip().split("|")[3])
            print each.strip().split("\t")[0]+"\t"+";".join(set(gene))#+"\t"+str(peplength)+"\t"+str(genespec)+"\t"+str(isoformspec)+"\t"+str(shared )+"\t"+";".join(nacc)
        
        
def count_isoforms_pergene(refseq):
    gene_isoform={}
    allgenes=[]
    #print "Gene\tCountofisoforms\tIsoforms"
    for each in SeqIO.parse(open(refseq),"fasta"):
        try:
            if each.id.split("|")[4].split("$")[0] in gene_isoform:
                gene_isoform[each.id.split("|")[4].split("$")[0]].append(each.id.split("|")[3])
            else:
                gene_isoform[each.id.split("|")[4].split("$")[0]]=[]
                gene_isoform[each.id.split("|")[4].split("$")[0]].append(each.id.split("|")[3])
        except:pass
    for gene,isoform in gene_isoform.iteritems():
        if len(isoform)>1:
            allgenes.append(gene)
        #print gene+"\t"+str(len(isoform))
    return allgenes
        
#genespecific_isoformspecific()



def map_peptides_to_refseq(fname):
    s=""
    acc=[]
    for i in SeqIO.parse(open("/home/srikanth/Downloads/CellMap/PSMs/Human_RefSeq65_withGS_061214.fasta"),"fasta"):
   	s+=str(i.seq)+"#"
   	acc.append(i.id)
    

    c=0
    nf=open("peptides_not_dbmatch.fasta","a")
    out=open("Peptides_matching_db"+fname.replace("./","")+".txt","w")
    for record in open(fname):
   	pep=record.strip()
   	match=[i.start() for i in re.finditer(pep.strip(),s)]
   	ind=[(acc[s[:j].count("#")]) for j in match]
   	if ind==[]:
  		nf.write(pep.strip()+"\t"+fname+"\n")
   	else:
  		out.write(pep.strip()+"\t"+";".join(ind)+"\n")
   	c+=1
   	if c%50==0:print "done",c
   	
map_peptides_to_refseq(sys.argv[1])