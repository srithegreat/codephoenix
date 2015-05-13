# -*- coding: utf-8 -*-
"""
Created on Mon Jan  5 18:32:50 2015
Modified 29th April
@author: srikanth
"""

#import pylab as P

import re,os,random,sys
from Bio import SeqIO
genome="/home/srikanth/Downloads/lncRNA_monoctyes/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
tcons_category={}
userinp=sys.argv[1].replace("./","")

  
def main(transcript):
    print "Working with ",userinp
    directory="temp"
    if not os.path.exists(directory):
        os.makedirs(directory)
    #os.system("cuffcompare -r RefSeq_Gencode_UCSC_cabili.gtf -o compared_ensembl_refseq_ucsc %s"%(transcript))
    #filter_gtf("compared_ensembl_refseq_ucsc."+userinp+".tmap","compared_ensembl_refseq_ucsc.combined.gtf")
       
    print "Creating fasta"
    temp="temp/temp_"+str(random.randint(1,100))+".fasta"
    #os.system("gffread %s -U -g %s -w %s"%("filtered_gtffile.gtf",genome,temp))#no single exon -U
    os.system("gffread %s -U -g %s -w %s"%(transcript,genome,temp))#Even single exon
    print "Done fasta, Filtering now...."
    filter_fasta(temp)
    print 'Done........'
    

def filter_gtf(tmap,gtf):
    
    #qual=filter_tmap(tmap)#tmap file
    fpkm=0.5
    qualified={}
    print "filtering class codes of interest"
    for i in open(tmap):
        if not i.startswith("ref_gene_id"):
            ref_gene_id,ref_id,class_code,cuff_gene_id,cuff_id,FMI,FPKM,FPKM_conf_lo,FPKM_conf_hi,cov,len,major_iso_id,ref_match_len=i.strip().split("\t")
            if (class_code=="u") and float(FPKM)>=fpkm:#class_code=="i" for lncRNA as against lincRNA
               # qualified.append(cuff_id)
               qualified[cuff_id]=FPKM+"#"+class_code
               
            
    filtered=open("temp/filtered_gtffile.gtf","w")
    for lines in open(gtf):
        cuffid=re.search("oId\s+\"(.*?)\"",lines.strip().split("\t")[8]).group(1)
        tcons=re.search("transcript_id\s+\"(.*?)\"",lines.strip().split("\t")[8]).group(1)
        if cuffid in qualified:
            tcons_category[tcons]=qualified[cuffid]
        if cuffid in qualified:
            filtered.write(lines)
            
    filtered.close()
    print "Written filtered gtf"
    
    

def parse_gtf(tcons):
    alllines=[]
    for lines in open("/home/srikanth/Downloads/lncRNA_monoctyes/cuffnorm_15samples/intergenic_intronic_transcripts_filtered_3x.gtf"):
       # cuffid=re.search("oId\s+\"(.*?)\"",lines.strip().split("\t")[8]).group(1)
        transcript_id=re.search("transcript_id\s+\"(.*?)\"",lines.strip().split("\t")[8]).group(1)
        if transcript_id in tcons:
           alllines.append(lines.strip()) 
    return alllines
    


            
    

def filter_fasta(fname):
    filtered_fasta=open("temp/filtered_transcripts.fasta","w")
    length=200
    records=0
    tc_seq={}
    for each in SeqIO.parse(open(fname),"fasta"):
        if len(each.seq)>length:
            records+=1
            filtered_fasta.write(">"+each.id+"\n"+str(each.seq)+"\n")
            tc_seq[each.id]=str(each.seq)
                   
            
    filtered_fasta.close()
    print "written records",records
    cpat_calc(tc_seq)
    
    
def cpat_calc(mapped):
    import numpy as np
    
    #get fpkm values from the isoforms_fpkm.table
    isoform_fpkm={}
    for iso in open("/home/srikanth/Downloads/lncRNA_monoctyes/cuffnorm_15samples/isoforms.fpkm_table"):
        if not iso.startswith("tracking"):
            isoform_fpkm[iso.strip().split("\t")[0]]="|".join(iso.strip().split("\t")[1:])+"|Avg_"+str(np.average([float(p) for p in iso.strip().split("\t")[1:16]]))
    print "calculating CPAT and filtering..."
    total=0
    out1=open("final_list_sorted"+userinp[:-4]+".txt","w")
    os.system("cpat.py -g temp/filtered_transcripts.fasta -o %s -x /media/srikanth/My\ Passport/backup_system/Downloads/CPAT-1.2.2/Human_Hexamer.tsv -d /media/srikanth/My\ Passport/backup_system/Downloads/CPAT-1.2.2/Human_logitModel.RData "%("temp/tempout"))
    filterid=[]
    #prob_transcripts_single=open("Putative_lncRNAs_Single_"+userinp[:-4]+".fasta",'w')
    prob_transcripts_multiple=open("Putative_lncRNAs_Multiple_"+userinp[:-4]+".fasta",'w')
    gtf="temp/"+userinp[:-4]+".gtf"
    
    for each in open("temp/tempout"):
        try:
            if float(each.strip().split("\t")[5])<0.375:
                total+=1
                filterid.append(each.strip().split("\t")[0])
        except:pass
    l=parse_gtf(filterid)
    op=open(gtf,'w')
    for j in l:
        op.write(j+"\n")
    op.close()
    #len_tcons=transcript_lengths(gtf)
    for i in filterid:
        prob_transcripts_multiple.write(">"+i+"|"+isoform_fpkm[i]+"\n"+mapped[i]+"\n")
        out1.write(i+"\t"+isoform_fpkm[i].split("Avg_")[1]+"\t"+mapped[i]+"\n")
    #print len_tcons
   
    #for i in filterid:
        
        #if len_tcons[i].startswith("Single"):
         #   prob_transcripts_single.write(">"+i+"|"+len_tcons[i]+"|"+tcons_category[i].split("#")[0]+"_Cat_"+tcons_category[i].split("#")[1]+"\n"+mapped[i]+"\n")
        #else:
        # out1.write(i+"\t"+len_tcons[i]+"\t"+tcons_category[i].split("#")[0]+"\t"+"Cat_"+tcons_category[i].split("#")[1]+"\t"+mapped[i]+"\n")
         #prob_transcripts_multiple.write(">"+i+"|"+len_tcons[i]+"|"+tcons_category[i].split("#")[0]+"_Cat_"+tcons_category[i].split("#")[1]+"\n"+mapped[i]+"\n")
    print "Total putative lncRNAs ",total
    #prob_transcripts_single.close()
    prob_transcripts_multiple.close()    
    
    



def transcript_lengths(gtf):
    transcripts={}
    multi=0
    single=0
    single_multi={}
    for lines in open(gtf):
        
        if lines.strip().split("\t")[2]=="exon":
            #print lines
            tcons=re.search("transcript_id\s+\"(.*?)\"",lines.strip().split("\t")[8]).group(1)
            if tcons+"$"+lines.strip().split("\t")[0] in transcripts:
                transcripts[tcons+"$"+lines.strip().split("\t")[0]].append((int(lines.strip().split("\t")[3]),int(lines.strip().split("\t")[4])+1))
                               
            else:
                transcripts[tcons+"$"+lines.strip().split("\t")[0]]=[]
                transcripts[tcons+"$"+lines.strip().split("\t")[0]].append((int(lines.strip().split("\t")[3]),int(lines.strip().split("\t")[4])+1))
                
                
                   
    for k,v in transcripts.iteritems():
        if len(v)>1:
            multi+=1
            single_multi[k.split("$")[0]]="Multi_"+str(sum([len(range(f[0],f[1]+1)) for f in v]))+"|"+k.split("$")[1]
            #print k
        else:
            single_multi[k.split("$")[0]]="Single_"+str(sum([len(range(f[0],f[1]+1)) for f in v]))+"|"+k.split("$")[1]
    #print multi
    
            
            
        
    #P.hist(lengths,bins=100)
    #P.pie([multi,len(transcripts)-single],labels=["Multi: "+str(multi),"Single: "+str(len(transcripts)-multi)])
    
    #P.show()
    #minlength=min(lengths)
    #maxlength=max(lengths)
    
    print "Single Exon: %s Multi-Exon: %s"%(len(transcripts)-multi,multi)
    #print "Min", minlength, "Max", maxlength
    #os.system("rm -rf compared_ensembl_refseq_ucsc.transcripts.gtf.tmap temp compared_ensembl_refseq_ucsc.combined.gtf")

    return single_multi


main(userinp)
