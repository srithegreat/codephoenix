##Program to map peptides to RefSeq or any protein fasta database##
import re
from Bio import SeqIO
s=""
acc=[]

for i in SeqIO.parse(open(raw_input("Enter fasta filename: ").strip()),"fasta"):
	s+=str(i.seq)+"#"
	acc.append(i.id)


c=0
nf=open("Peptides_not_mapping_refseq.fasta","w")#output filename with unique entries
out=open("Peptides_matching_all.txt","w") #output filename with_all mappings
#pepcol=int(raw_input("Enter Column of peptides"))-1
for record in open(raw_input("Enter filename with peptides: ").strip()):
	
	pep=record.strip()
	#header=record.description
	match=[i.start() for i in re.finditer(pep.strip().upper(),s)]
	#print match
	ind=[(acc[s[:j].count("#")]) for j in match]
	#print pep,ind
	#if ind==[]:
	#nf.write(">"+header+"\n"+pep+"\n")
	#else:
	out.write(pep.strip()+"\t"+";".join(ind)+"\n")
	c+=1
	if c%50==0:print "done",c


