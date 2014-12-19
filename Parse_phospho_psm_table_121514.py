# -*- coding: utf-8 -*-
"""
Created on Mon Dec 15 13:16:53 2014

@author: srikanth

Program to parse phospho experiments

Input text file with PSMs:
Sequence	Modifications	phosphoRS Site Probabilities	Spectrum File
"""
import re,sys
from Bio import SeqIO
pepglobal={}
import argparse

def correct_assignment_Y(database,inputfile):
   """Program only corrects sites for pTyr"""
   psm_site={}
   line=0
   out=open(inputfile[:-4]+"_corrected.txt",'w')
   out1=open("PSMs_Corrected_mapped_to_proteins.txt","w")
   #out.write("Status\t#of Phosphorylation Events detected	Peptide	Corrected Site (75%)	Sequence	Modification (s)	Probability	\n")
   cutoff=75.0
   
   for eachline in open(inputfile):
       
       if not eachline.startswith("Sequence"):
           line+=1
           
           #if len(eachline.split("\t"))==4 and eachline.split("\t")[2]!="":
           peptide,modification,phospho=eachline.strip().split("\t")[0],eachline.strip().split("\t")[1],eachline.strip().split("\t")[2]#specify columns from the file
           if "Too" in phospho:
               continue
           if phospho=="":
               continue
           original=peptide.replace("m","M").replace("a","A").replace("c","C").replace("r","R").replace("k","K").replace("k","K").replace("g","G").replace("v","V")
           allphosmod=re.findall("([STY]\d+)\(Phospho\)",modification)
           if allphosmod!=[]:
               sequence=peptide.upper()
               phosphoRS=[s.replace("(","").replace(")","") for s in phospho.split(";")]
               f={}
               noofphospho=len(allphosmod)
               for e in phosphoRS:
                   f[e.split(":")[0].strip()]=float(e.split(":")[1].strip())
               qualifiedsites=[]
               alldict=[]
               
               for p,v in f.iteritems():
                   if float(v)==max(f.values()):
                       alldict.append(p)
                   if float(v)>=cutoff:
                       #print p,v,sequence,allphosmod
                       qualifiedsites.append(p)
               newsequence=[i.upper() for i in list(sequence)]
               for every in qualifiedsites:
                   psm_site[peptide.upper()+"_"+"#".join(qualifiedsites)]=0
                   residue=re.match("(\w)(\d+)",every).group(1)
                   pos=re.match("(\w)(\d+)",every).group(2)
                   #newsequence=sequence[:int(every[1])-1].replace("s","S").replace("t","T").replace("y","Y")+every[0].lower()+sequence[int(every[1]):].replace("s","S").replace("t","T").replace("y","Y")
                   newsequence[int(pos)-1]=residue.lower()
               status="Corrected"
               if qualifiedsites!=[]:
                   if original=="".join(newsequence):
                       status="Not Corrected"
                       
                   out.write("\t".join([status+"\t"+str(noofphospho),"".join(newsequence),";".join(qualifiedsites)])+"\t"+eachline.strip()+"\n")
               else:
                   status="Corrected"
                   corsite=[ns for ns in alldict]
                   corsiteY=[ns for ns in alldict if ns.startswith("Y")]
                   
                   
                   if len(corsite)>1:
                       corsite=allphosmod
                       status="Now Corrected"
                       for every2 in corsite:
                           psm_site[peptide.upper()+"_"+"#".join(corsite)]=0
                           residue=re.match("(\w)(\d+)",every2).group(1)
                           pos=re.match("(\w)(\d+)",every2).group(2)
                           #newsequence=sequence[:int(every[1])-1].replace("s","S").replace("t","T").replace("y","Y")+every[0].lower()+sequence[int(every[1]):].replace("s","S").replace("t","T").replace("y",Y")
                           
                           newsequence[int(pos)-1]=residue.lower()
                   if corsite==corsiteY and len(corsite)>1:
                       status="Ambiguous Y"
                       print corsiteY,corsite
                   if len(corsiteY)==0:
                       status="Ambiguous ST"
                       newsequence=sequence
                   out.write("\t".join([status+"\t"+str(noofphospho),"".join(newsequence),";".join(corsite)])+"\t"+eachline.strip()+"\n")
                   
   out.close()
   mappedpeps=map_proteins(psm_site,database)
    #print mappedpeps.keys()
   
   for j in open(inputfile[:-4]+"_corrected.txt"):
       
       #print i.strip().split("\t")[1]
  
       if mappedpeps.has_key(j.strip().split("\t")[2].upper()+"_"+j.strip().split("\t")[3]):
         out1.write(j.strip()+"\t"+mappedpeps[j.strip().split("\t")[2].upper()+"_"+j.strip().split("\t")[3]]+"\n")
       else:
          out1.write(j.strip()+"\t"+"-\n")
       
   out1.close()
 


def correct_assignment_ST():
   """Program only corrects sites for pST"""
   psm_site={}
   line=0
   out=open(sys.argv[1][:-4]+"_corrected.txt",'w')
   out1=open("PSMs_Corrected_mapped_to_proteins.txt","w")
   out.write("Status\t#of Phosphorylation Events detected	Peptide	Corrected Site (75%)	Sequence	Modification (s)	Probability\n")
   cutoff=75.0
   
   for eachline in open(sys.argv[1]):
       
       if not eachline.startswith("Sequence"):
           line+=1
           #print eachline
           
           peptide,modification,phospho=eachline.strip().split("\t")[0],eachline.strip().split("\t")[1],eachline.strip().split("\t")[2]#specify columns from the file
           if "Too" in phospho:
               continue
           original=peptide.replace("m","M").replace("a","A").replace("c","C").replace("r","R").replace("k","K").replace("k","K").replace("g","G").replace("v","V")
           allphosmod=re.findall("([STY]\d+)\(Phospho\)",modification)
           if allphosmod!=[]:
               sequence=peptide.upper()
               phosphoRS=[s.replace("(","").replace(")","") for s in phospho.split(";")]
               f={}
               noofphospho=len(allphosmod)
               for e in phosphoRS:
                   f[e.split(":")[0].strip()]=float(e.split(":")[1].strip())
               qualifiedsites=[]
               alldict=[]
               
               for p,v in f.iteritems():
                   if float(v)==max(f.values()):
                       alldict.append(p)
                   if float(v)>=cutoff:
                       #print p,v,sequence,allphosmod
                       qualifiedsites.append(p)
               newsequence=[i.upper() for i in list(sequence)]
               for every in qualifiedsites:
                   psm_site[peptide.upper()+"_"+"#".join(qualifiedsites)]=0
                   residue=re.match("(\w)(\d+)",every).group(1)
                   pos=re.match("(\w)(\d+)",every).group(2)
                   #newsequence=sequence[:int(every[1])-1].replace("s","S").replace("t","T").replace("y","Y")+every[0].lower()+sequence[int(every[1]):].replace("s","S").replace("t","T").replace("y","Y")
                   newsequence[int(pos)-1]=residue.lower()
               status="Corrected"
               if qualifiedsites!=[]:
                   if original=="".join(newsequence):
                       status="Not Corrected"
                       
                   out.write("\t".join([status+"\t"+str(noofphospho),"".join(newsequence),";".join(qualifiedsites)])+"\t"+eachline.strip()+"\n")
               else:
                   status="Ambiguous"
                   
                   out.write("\t".join([status+"\t"+str(noofphospho),"".join(original),";".join(allphosmod)])+"\t"+eachline.strip()+"\n")
   
   mappedpeps=map_proteins(psm_site)
    #print mappedpeps.keys()
   
   for j in open(sys.argv[1][:-4]+"_corrected.txt"):
       
       #print i.strip().split("\t")[1]
       if not j.startswith("#"):
           if mappedpeps.has_key(j.strip().split("\t")[2].upper()+"_"+j.strip().split("\t")[3]):
             out1.write(j.strip()+"\t"+mappedpeps[j.strip().split("\t")[2].upper()+"_"+j.strip().split("\t")[3]]+"\n")
           else:
              out1.write(j.strip()+"\t"+"-\n")
           
   out1.close()              

def map_proteins(peptides,database):
    print "working ",
    sys.stdout.flush()
    s=""
    acc=[]
    seqs={}
    #for i in SeqIO.parse(("/home/srikanth/Downloads/Refseq63_stats/human.protein.acc.RefSeq63.021014.fasta"),"fasta"):"/media/srikanth/srikanth_HDD1/IL-33/mouse_cont_RefSeq60_JUL2013.fasta"
    for i in SeqIO.parse(database,"fasta"):
        if not i.id.split("|")[1].startswith("c"):#remove contaminants
            s+=str(i.seq)+"#"
            acc.append(i.description)
            seqs[i.id.split("|")[1]]=str(i.seq)
    
    c=0
    out=open("peptides_matched_temp.txt","w")
    out.write("Peptide\tGene(s)\tPeptidePosition(s)\tProteinAccession(s)\tProteinLength(s)\n")
    for p in peptides:
        pep=p.split("_")[0]
    
        match=[i.start() for i in re.finditer(pep.strip().upper(),s)]
        
        #print match
        
        ind=[(acc[s[:j].count("#")]) for j in match]
        #print pep,ind
        #print pep,ind
        ids=[l.split("|")[1] for l in ind]
        #protpos=[str(seqs[k].index(pep)+1) for k in ids]
        protlengths=[str(len(seqs[k])) for k in ids]
        protacc=[l.split("|")[3] for l in ind]
        peppos=""
        for pos in p.split("_")[1].split("#"):
            peppos+=";".join([pos+"_"+str(seqs[k].index(pep)+int(re.match("\w(\d+)",pos).group(1))) for k in ids])+"|"
        pepglobal[p]=";".join([l.split("|")[4].split("$")[1] for l in ind])+"\t"+peppos[:-1]+"\t"+";".join(protacc)+"\t"+";".join(protlengths)
        out.write(pep.strip()+"\t"+";".join([l.split("|")[4].split("$")[0] for l in ind])+"\t"+peppos[:-1]+"\t"+";".join(protacc)+"\t"+";".join(protlengths)+"\n")
        c+=1
        
        if c%50==0:
            print "\b*",
            sys.stdout.flush()
    
    out.close()
    return pepglobal   
    
def create_sitelist():
    nf=open("Peptide_site_allidentified_othergenes.txt","w")
    nf.write("Peptide	AllGenes	AllAccessions	Leadingprotein(length)	Accession	Gene Symbol	ProteinPosition	Status	Other Sites	Othergenes\n")
    fprint=set([])
   
    for each in open("peptides_matched_temp.txt"):
        if not each.startswith("Pep"):
            lines=each.strip().split("\t")
            if len(lines)>1:
                gene,proteinpos,protein,length=lines[1],lines[2],lines[3],lines[4]
                lenlist=[int(k) for k in length.split(";")]
                indmax=lenlist.index(max([int(p) for p in lenlist]))
                
                genelist=gene.split(";")
                proteinlist=protein.split(";")
                nglist=list(set([s for s in genelist if s!=genelist[indmax]]))
                allnps=[f3 for f3 in proteinlist if f3.startswith("NP")]
                allxps=[f4 for f4 in proteinlist if f4.startswith("XP")]
                lengthnps=[lenlist[proteinlist.index(f5)] for f5 in allnps]
                
                longest=proteinlist[indmax]
                if longest.startswith("NP"):
                    longest=longest
                elif longest.startswith("XP") and allnps==[]:
                    longest=longest   
                    
                elif longest.startswith("XP") and allnps!=[]:
                    longest=allnps[lengthnps.index(max(lengthnps))]
                    indmax=proteinlist.index(allnps[lengthnps.index(max(lengthnps))])               
                    
                
                    
                
                
                if len(set(genelist))!=1:
                    status="Shared"
                    nglist=list(set([s for s in genelist if s!=genelist[indmax]]))
                    posforother=[proteinpos.split(";")[genelist.index(k)] for k in nglist]
                    if proteinpos.count("|")>0:
                        
                        
                        prottemp=""
                        for i1 in proteinpos.split("|"):
                            prottemp+=re.search("([STY])\d+_(\d+)",i1.split(";")[indmax]).group(1)+re.search("([STY])\d+_(\d+)",i1.split(";")[indmax]).group(2)+";"
                        #fprint.add("\t".join([lines[0],lines[1],lines[3],str(lenlist[indmax]),proteinlist[indmax],genelist[indmax],prottemp,status,";".join(posforother),";".join(nglist)]))
                        fprint.add("\t".join([lines[0],lines[1],lines[3],str(lenlist[indmax]),"".join(longest),genelist[indmax],prottemp,status,";".join(posforother),";".join(nglist)]))
                    else:
                        proteinposlist=proteinpos.split(";")
                        
                        fprint.add("\t".join([lines[0],lines[1],lines[3],str(lenlist[indmax]),"".join(longest),genelist[indmax],re.search("([STY])\d+_(\d+)",proteinposlist[indmax]).group(1)+re.search("([STY])\d+_(\d+)",proteinposlist[indmax]).group(2),status,";".join(posforother),";".join(nglist)]))
                    #print longest
                    
                else:
                    status="Unique"
                    if proteinpos.count("|")>0:
                        prottemp=""
    
                        for i1 in proteinpos.split("|"):
                            prottemp+=re.search("([STY])\d+_(\d+)",i1.split(";")[indmax]).group(1)+re.search("([STY])\d+_(\d+)",i1.split(";")[indmax]).group(2)+";"
                        #fprint.add("\t".join([lines[0],lines[1],lines[3],str(lenlist[indmax]),proteinlist[indmax],genelist[indmax],prottemp,status]))
                        fprint.add("\t".join([lines[0],lines[1],lines[3],str(lenlist[indmax]),"".join(longest),genelist[indmax],prottemp,status]))
                    else:
                        proteinposlist=proteinpos.split(";")
                        
                        fprint.add("\t".join([lines[0],lines[1],lines[3],str(lenlist[indmax]),"".join(longest),genelist[indmax],re.search("([STY])\d+_(\d+)",proteinposlist[indmax]).group(1)+re.search("([STY])\d+_(\d+)",proteinposlist[indmax]).group(2),status]))

    for i in list(fprint):
       nf.write(i+"\n")
            

def getgenewise():
    fo=open("Peptide_gene_Protsite_nr.txt",'w')  
    uniqdict={}          
    for each in open("peptides_matched_temp.txt"):
        try:
            if not each.startswith("Peptide"):
                peptide,gene,proteinpos=each.strip().split("\t")[0],each.strip().split("\t")[1],each.strip().split("\t")[2]
                cnt=0
                for i in gene.split(";"):
                    if proteinpos.count("|")==0:
                        uniqdict["\t".join([peptide,i,re.search("([STY])\d+_(\d+)",proteinpos.split(";")[cnt]).group(1),re.search("([STY])\d+_(\d+)",proteinpos.split(";")[cnt]).group(2)])]=0
                        #fo.write("\t".join([peptide,i,re.search("([STY])\d+_(\d+)",proteinpos.split(";")[cnt]).group(1),re.search("([STY])\d+_(\d+)",proteinpos.split(";")[cnt]).group(2)])+"\n")
                        cnt+=1
                    elif proteinpos.count("|")>0:
                        for k1 in proteinpos.split("|"):
                            uniqdict["\t".join([peptide,i,re.search("([STY])\d+_(\d+)",k1.split(";")[cnt]).group(1),re.search("([STY])\d+_(\d+)",k1.split(";")[cnt]).group(2)])]=0
                            #fo.write("\t".join([peptide,i,re.search("([STY])\d+_(\d+)",k1.split(";")[cnt]).group(1),re.search("([STY])\d+_(\d+)",k1.split(";")[cnt]).group(2)])+"\n")
                        cnt+=1
        except:pass
    for key in uniqdict:
        fo.write(key+"\n")
def uniq():
    """Generated unique Gene and all sites Use after getgenewise()"""
    nf=open("Genes_sites.txt","w")
    nf.write("Gene Symbol\tSite(s)\n")
    f={}
    for i in open("Peptide_gene_Protsite_nr.txt"):
        if f.has_key(i.strip().split("\t")[1]):
            f[i.strip().split("\t")[1]].append(i.strip().split("\t")[2]+i.strip().split("\t")[3])
        else:
            f[i.strip().split("\t")[1]]=[]
            f[i.strip().split("\t")[1]].append(i.strip().split("\t")[2]+i.strip().split("\t")[3])
    
    for k,v in f.iteritems():
        temp=[]
        for l in v:
            if not l in temp:
                temp.append(l)
        
        nf.write( k+"\t"+";".join(temp)+"\n")




def getpsmcount_labelfree():
    nf1=open("Final_peptides_with_PSM_counts.txt","w")
    f4={}
    for each in open("peptides_matched_temp.txt"):
        if not each.startswith("Peptide"):
           #print each
            peptide,gs,site,prot,length=each.strip().split("\t")
            if site.count("|")==0:
                
                f4[peptide+"#"+re.search("([STY]\d+)_(\d+)",site).group(1)]=";".join(set(gs.split(";")))+"\t"+";".join(set(prot.split(";")))
            else:
                newsite=""
                for every in site.split("|"):
                    newsite+=re.search("([STY]\d+)_(\d+)",every).group(1)+";"
                
                f4[peptide+"#"+newsite[:-1]]=gs+"\t"+prot
		
		 

    #for each in open("temp"):
     #   if f.has_key(each.strip().split("\t")[0]+"#"+each.strip().split("\t")[1]):
      #      print f[each.strip().split("\t")[0]+"#"+each.strip().split("\t")[1]]+"\t"+each.strip()    
    
    
    g=[]
    for f1 in open("PSMs_Corrected_mapped_to_proteins.txt"):
        if not f1.startswith("Status"):
            
            if not f1.strip().split("\t")[7] in g and f1.strip().split("\t")[7].endswith("raw"):
                g.append(f1.strip().split("\t")[7])
                

    f={}
    for each in open("PSMs_Corrected_mapped_to_proteins.txt"):
        if not each.startswith("Status"):
            peptide,site,raw=each.strip().split("\t")[2],each.strip().split("\t")[3],each.strip().split("\t")[7]
            if f.has_key(peptide.upper()+"$"+site):
                f[peptide.upper()+"$"+site].append(raw)
            else:
                f[peptide.upper()+"$"+site]=[]
                f[peptide.upper()+"$"+site].append(raw)
          
    nf1.write( "GeneSymbol (s)\tAccession(s)\tPeptide\tSite\t"+"\t".join(g)+"\n")
    for k,v in f.iteritems():
        r={}
        
        for raw in v:
            if r.has_key(raw):
                r[raw]+=1
            else:
                r[raw]=1
        out=""
        for l in g:
            if r.has_key(l):
                out+=str(r[l])+"\t"
            else:
                out+="0\t"
                
        #print k.split("$")[0]+"\t"+k.split("$")[1]+"\t"+out
        nf1.write( f4[k.split("$")[0]+"#"+k.split("$")[1]]+"\t"+k.split("$")[0]+"\t"+k.split("$")[1]+"\t"+out+"\n")
        
#uniq()

def writexlsx():
    """Writes final file in xlsx format"""
    import xlsxwriter
    newxls=xlsxwriter.Workbook("Genes_sites_Psmcount.xlsx")
    worksheet=newxls.add_worksheet()
    row=0
    #col=0
    for l in open("Final_peptides_with_PSM_counts.txt"):
        lenlines=len(l.strip().split("\t"))
        for i in range(lenlines):
            worksheet.write(row,i,l.strip().split("\t")[i])
        row+=1
    newxls.close()
    
if __name__=="__main__":
    
    __author__ = 'Srikanth'
    parser = argparse.ArgumentParser(description='Program to parse PSM list, correct phospho sites and print PSM counts.')
    parser.add_argument('-d','--database', help='Input database name',required=True)
    parser.add_argument('-i','--inputfile', help='Input file',required=True)
    #parser.add_argument('-o','--output',help='Output file name', required=True)
    args = parser.parse_args()

    correct_assignment_Y(args.database,args.inputfile)
    create_sitelist()
    getgenewise()
    uniq()
    getpsmcount_labelfree()
    writexlsx()