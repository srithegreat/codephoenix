# -*- coding: utf-8 -*-
"""
Created on Fri Oct 18 17:12:04 2013

@author: srikanth
"""

import argparse
import glob
import numpy as np
from math import sqrt




def parse(inpfile,tum):
    if tum==True:
        tempoutput="T_temp.txt"
    elif tum==False:
        tempoutput="N_temp.txt" 
    alltxt=glob.glob(inpfile+"/*.txt")#Full path of files
    output=open(tempoutput,'w')#outputfile name
    probescpg={}
    probe_value={}
    for p in open("/home/srikanth/Downloads/TCGA/probes_in_TSS_Island_region.tsv"):#probes_in_TSS_only_region.tsv"):
    	probescpg[p.strip().split("\t")[0]]=p.strip().split("\t")[21]
    
    	
    count=1
    for f in alltxt:
        print "reading file", count
        count+=1
        c=0
        for each in open(f):
            if not each.startswith("probe"):
                probe,betavalue=each.strip().split("\t")[0],each.strip().split("\t")[1]#check the exact column
                if betavalue!="NA":
                    c+=1
                    if probescpg.has_key(probe):
                        if probe_value.has_key(probe):
                            probe_value[probe].append(float(betavalue))
                        else:
                            probe_value[probe]=[]
                            probe_value[probe].append(float(betavalue))
    	#print len(probe_value)						
    	cnt=0
    for p,beta in probe_value.iteritems():
        cnt+=1
        #print p+"\t"+"\t".join([str(b) for b in beta])+"\n"
        output.write(p+"\t"+";".join(set(probescpg[p].split(";")))+"\t"+"\t".join([str(b) for b in beta])+"\n")         
    print "calculated sample mean and SD..grouping gene wise probes.."
    output.close()
    group_genewise_probes(tempoutput,tum)
      
def group_genewise_probes(temp,tum):
    """Averages probes in TSS region genewise"""
    if tum==True:
        output="T_tempG.txt"
        f=0
    else:
        output="N_tempG.txt"
        f=1
    out=open(output,'w')
    c=0
    
    genebeta={}
    for each in open(temp):#file containes probes and beta values
    	probe,genes,betavalues=each.strip().split("\t")[0],each.strip().split("\t")[1],[float(s) for s in each.strip().split("\t")[2:]]
    	
    	if len(genes.split(";"))==1:
    		#print genes
    		if genebeta.has_key(genes):
            		genebeta[genes].append(betavalues)
    		else:
    			genebeta[genes]=[]
    			genebeta[genes].append(betavalues)
    #print genebeta
    for gene,probes in genebeta.iteritems():
            c+=1
            s_mean=map(np.mean,zip(*probes))
            gene_mean,gene_std=np.mean(s_mean),np.std(s_mean)
            out.write( gene+"\t"+str(gene_mean)+"\t"+str(gene_std    )+"\n")
    print "Done genewise averaging tumor.."
    print "working with normals.."
    
    if f==0:
        parse(args.locationnormal,tum=False)
        
        print "Done genewise averaging normal.."
        print "Comparing tumor normal"
        compare_tum_norm()
        f=1
    

def compare_tum_norm():
   """ if tumor mean>2sd normal mean"""
   k={}
   output=args.output
   nfile=open(output,"w")
   nfile.write("Genesymbol\tTumor_mean\tTumor_SD\tNormal_mean\tNormal_SD\tStatus\n")
   for i in open("T_tempG.txt"):#Tumor
       
       
       k[i.strip().split("\t")[0]]=(i.strip().split("\t")[1])+"\t"+(i.strip().split("\t")[2])
   
   for j in open("N_tempG.txt"):#Normal
       status="NA"
       if k.has_key(j.strip().split("\t")[0]):
           if float(k[j.strip().split("\t")[0]].split("\t")[0])>float(j.strip().split("\t")[1])+2*float(j.strip().split("\t")[2]):
               status="hyper"
               nfile.write( j.strip().split("\t")[0]+"\t"+k[j.strip().split("\t")[0]]+"\t"+j.strip().split("\t")[1]+"\t"+j.strip().split("\t")[2]+"\t"+status+"\n")
           elif float(j.strip().split("\t")[1])>float(k[j.strip().split("\t")[0]].split("\t")[0]):
               status="hypo"
               nfile.write( j.strip().split("\t")[0]+"\t"+k[j.strip().split("\t")[0]]+"\t"+j.strip().split("\t")[1]+"\t"+j.strip().split("\t")[2]+"\t"+status+"\n")

       #else:
        #   nfile.write( j.strip().split("\t")[0]+"\t"+"NA\tNA"+"\t"+j.strip().split("\t")[1]+"\t"+j.strip().split("\t")[2]+"\t"+status+"\n")
   print "Done..Thanks you!!!!"
def meanstdv(x):
	
	n, mean, std = len(x), 0, 0
	for a in x:
		mean = mean + a
	mean = mean / float(n)
	for a in x:
		std = std + (a - mean)**2
	std = sqrt(std / float(n))#changed SD to PD
	return mean, std

	
if __name__=="__main__":
    
    __author__ = 'Srikanth'
     
    parser = argparse.ArgumentParser(description='Program to parse TCGA methylation data and create hypermethylated gene list.')
    parser.add_argument('-t','--locationtumor', help='Input file(s) location tumor',required=True)
    parser.add_argument('-n','--locationnormal', help='Input file(s) location normal',required=True)
    parser.add_argument('-o','--output',help='Output file name', required=True)
    args = parser.parse_args()

    parse(args.locationtumor,tum=True)
    
