#requirement
#bPTR (https://github.com/christophertbrown/iRep/blob/master/iRep/bPTR)
import sys,os

if len(sys.argv)<4:
        print("USAGE: python identify_PTR_from_pathoid.py sam_report_tsv shrunk_samfile_BEFORE_pathoid outDir")
else:
        sam_report_tsv=sys.argv[1]
        mappingfile="/mnt/bioinfonew/srikanth/Test_align/genomeindices/temp_genomes/mapping_genomes.txt"
        samfile=sys.argv[2]

        outdir=sys.argv[3]
        map_dict={}
        if not os.path.exists(outdir):
            os.makedirs(outdir)
         #os.system("mkdir %s"%outdir)
        for i in open(mappingfile):
                map_dict[i.strip().split("\t")[1]]=i.strip().split("\t")[0]

        for k in open(sam_report_tsv):
                if not (k.startswith("Total") or k.startswith("Genome")) :
                        mapped_genomeid=k.strip().split("\t")[0]
                        fileid=map_dict[mapped_genomeid]
                        print("Running PTR prediction for",mapped_genomeid)
                        os.system("bPTR -f /mnt/bioinfonew/srikanth/Test_align/genomeindices/temp_genomes/%s -t 8 -s %s -o %s/bPTR_%s.tsv -plot %s/bPTR_%s.pdf -m coverage"%(fileid,samfile,outdir,fileid,outdir,fileid))
