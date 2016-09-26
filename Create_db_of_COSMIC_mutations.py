##Program to create mutant protein db from COSMIC db###


from Bio import SeqIO
import re
mutations={}
dictofnames={}
def getinfo():
	for each in open("COSMIC66_MART.txt"):#Replace with latest cosmic mart file
		try:
			##	print each
			dictofnames[each.strip().split("\t")[4]+"#"+each.strip().split("\t")[16].split(".")[1]]=each.strip().split("\t")[15]
		except:pass

	
		
def parse_file():
	doubt=open("ambiguous_cases_nonsynonymous.txt",'w')
	counttryp=0
	mutcount=0
	coordinatelist={}
	for mutation in open("non_synonymous_mutations_COSMIC.txt"):##non_synonymous_mutations_COSMIC.txt
		try:
			gene,amino_residue=mutation.strip().split("\t")[0],mutation.strip().split("\t")[14].replace("*","Z")
			aminowild=re.search("p.(\w)(\d+)(\w)",amino_residue).group(1)
			aminomutant=re.search("p.(\w)(\d+)(\w)",amino_residue).group(3)
			position=re.search("p.(\w)(\d+)(\w)",amino_residue).group(2)
			if aminowild=="Z":aminowild="*"
			if aminomutant=="Z":aminomutant="*"
			
			if coordinatelist.has_key(gene):
				coordinatelist[gene][position]=(aminowild,aminomutant)
				#coordinatelst[gene][position].append(aminowild)
			else:
				coordinatelist[gene]={}
				coordinatelist[gene][position]=(aminowild,aminomutant)
				#coordinatelst[gene][position].append(aminowild)
	
	
		except:
			
			doubt.write(mutation.strip()+"\n")

	replaceamino(coordinatelist)

def replaceamino(coordinates):
	#regex = re.compile("([KR]?[^P].*?[KR](?!P))")
		
	mutants=open("mutated_sequences.fasta","w")
	mutant_tryptic=open("mutant_sequences_tryptic.fasta","w")
	n=0
	for record in SeqIO.parse(open("allproteins_cosmic.fasta"),"fasta"):
		sequence=str(record.seq).replace("*","")
		
		if coordinates.has_key(record.id):
			positions=coordinates[record.id]
			for pos,change in positions.iteritems():
				try:
					pos_in_genome=dictofnames[record.id+"#"+change[0]+str(pos)+change[1]]
				except:pass
				mutantsequence=sequence[:int(pos)-1]+change[1]+sequence[int(pos):]
				temp=re.sub("(?<=[KR])(?!P)","\n",mutantsequence)
				tryp=temp.split()
				
				for each in tryp:
					i=mutantsequence.index(each)+1
					c=0
					#print i,pos
					if (len(each)>5 and len(each)<26):
						
						if int(pos) in xrange(i,i+len(each)):
							n+=1
							
							peppos=int(pos)-(i)
							header=">COS"+str(n)+"M"+"|"+change[0]+str(pos)+change[1]+"|pep."+str(peppos)+"|"+record.id+"|"+pos_in_genome
							#final=[pep for pep in each.split("*") if len(pep)==max(len(each.split("*")))]
									
							mutant_tryptic.write(header+"\n"+each+"\n")
							
						
							
	
				
			
			






						
			
				

if __name__=="__main__":
	getinfo()	
	parse_file()
