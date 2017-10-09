from Bio import SeqIO
import sys
import time
start=time.clock()
f={}
c=0
for i in SeqIO.parse(open(sys.argv[1]),"fasta"):
	c+=1
	#print c
	if str(i.seq) in f:
		f[str(i.seq)].append(i.description)
	else:
		f[str(i.seq)]=[]
		f[str(i.seq)].append(i.description)
end=time.clock()
print "Done reading in ",(end-start)/60

fn=open("nrfile.fasta",'w')#file with non-redundant sequences
head=open("nrfile.fasta.index",'w') #file with index
count=0
for key,value in f.iteritems():
	count+=1
	header=">SFG"+str(count)
	head.write(header+"\t"+",".join(value)+"\n")
	fn.write(header+"\n"+key+"\n")





