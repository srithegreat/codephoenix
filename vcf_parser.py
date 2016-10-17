#requires pyvcf https://pyvcf.readthedocs.io/en/latest/#
import vcf
vcf_reader = vcf.Reader(open('IOB_0005.hg19_multianno_BLOOD.vcf', 'r'))
nonsynonymous=0
exonic=0
nonsynonymousq=0
exonicq=0
total=0
cat={}


for record in vcf_reader:
	total+=1
	
		
	
	#if record.INFO["ExonicFunc.refGene"][0]=="nonsynonymous_SNV":
	#	nonsynonymous+=1
	#elif record.INFO["Func.refGene"][0]=="exonic":
	#	exonic+=1
	for s in record.samples:
		if s["AD"][1]>=10:
			if record.INFO["Func.refGene"][0] in cat:
				cat[record.INFO["Func.refGene"][0]]+=1
			else:
				cat[record.INFO["Func.refGene"][0]]=0
				cat[record.INFO["Func.refGene"][0]]+=1
			if record.INFO["ExonicFunc.refGene"][0]=="nonsynonymous_SNV":
				nonsynonymousq+=1
			elif record.INFO["Func.refGene"][0]=="exonic":
				exonicq+=1
print "total", total 
print "qualifying", nonsynonymousq
print "Qualifying exonic", exonicq
for i,v in cat.iteritems():
	print i,v
