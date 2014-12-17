# -*- coding: utf-8 -*-
"""
Created on Fri Apr  4 19:41:43 2014
Getflanking amino acid residue for motif x
INPUT: txt file with peptide, accession, proteinposition
@author: srikanth
"""
import sys,re
from Bio import SeqIO
def getflanking():
    s={}
    for record in SeqIO.parse(open(sys.argv[1]),"fasta"):
        s[record.id.split("|")[3]]=str(record.seq)
        
        
    lineno=0
    for each in open(sys.argv[2]):
        if not each.startswith("#"):
            lineno+=1
            peptide,accession,position=each.strip().split("\t")[0],each.strip().split("\t")[1],each.strip().split("\t")[2]
            protpos=re.search("[STY](\d+)",position).group(1)
            residue=re.search("([STY])\d+",position).group(1)
            if position.count(";")==0:
                if s.has_key(accession):
                    sequence=s[accession]
                    up=sequence[int(protpos)-8:int(protpos)-1]
                    down=sequence[int(protpos):int(protpos)+7]
                    length=len(sequence)
                    if int(protpos)<8 :
                        
                        up=sequence[:int(protpos)-1]
                        numberofdash="_"*(7-len(up))
                        down=sequence[int(protpos):int(protpos)+7]
                        print numberofdash+up+residue.lower()+down+"\t"+each.strip()
                    elif int(protpos)>(length-8):
                        
                        up=sequence[int(protpos)-8:int(protpos)-1]
                        
                        down=sequence[int(protpos):]
                        numberofdash="_"*(7-len(down))
                        print up+residue.lower()+down+numberofdash+"\t"+each.strip()+"\t"+str(lineno)
                    else:
                        print up+residue.lower()+down+"\t"+each.strip()+"\t"+str(lineno)
            else:
                for l in position.split(";")[:-1]:
                    k=re.search("[STY](\d+)",l).group(1)
                    if s.has_key(accession):
                        sequence=s[accession]
                        up=sequence[int(k)-8:int(k)-1]
                        down=sequence[int(k):int(k)+7]
                        length=len(sequence)
                        if int(k)<8 :
                            
                            up=sequence[:int(k)-1]
                            numberofdash="_"*(7-len(up))
                            down=sequence[int(k):int(k)+7]
                            print numberofdash+up+residue.lower()+down+"\t"+each.strip()+"\t"+str(lineno)
                        elif int(k)>(length-8):
                            
                            up=sequence[int(k)-8:int(k)-1]
                            
                            down=sequence[int(protpos):]
                            numberofdash="_"*(7-len(down))
                            print up+residue.lower()+down+numberofdash+"\t"+each.strip()+"\t"+str(lineno)
                        else:
                            print up+residue.lower()+down+"\t"+each.strip()+"\t"+str(lineno)

if len(sys.argv)<2:
    print "USAGE: getflanking.py databasefasta txtfile"
else:                
    getflanking()
