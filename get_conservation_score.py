import sys, os, re
import numpy as np
if len(sys.argv)<2:
	print "Usage: python get_conservation_score.py gtffilename"
else:

	gtf=sys.argv[1]# gtf file
	nu=1
	for i in open(gtf):
		if not i.startswith("#"):
			lines=i.strip().split("\t")
			if lines[2]=="exon":
				transcriptid= re.search("transcript_id\s+\"(.*?)\"",lines[8]).group(1)
				chromosome, start, end= lines[0], int(lines[3]), int(lines[4])
				os.system("./bigWigToBedGraph hg38.phastCons20way.bw transcriptid.bedGraph -chrom=%s -start=%d -end=%d"%(chromosome, start, end))
				total=[]
				try:
					for l in open("transcriptid.bedGraph"):
						chromosome,start,end,score=l.strip().split("\t")
						total.append(float(score))
					print transcriptid, np.average(total)
				except:print transcriptid, "0"
