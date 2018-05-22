## Pre-requisites ##
#python3 or above
#Gemtools, GEM mapper (https://github.com/gemtools/gemtools)
#Trimgalore (https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)
#shrinksam (https://github.com/bcthomas/shrinksam)
#Pathoscope (https://github.com/PathoScope/PathoScope)

#usage: time sh run_pathomap_ID_pipeline.sh <GENOME.FASTA> <forward.fastq> <reverse.fastq> <Directoryname>

#check if directory exists
if [ ! -d "${4}" ]; then
mkdir ${4}
fi


#index new genome with name indexedgenome 
# gem-indexer -i {1} -o indexedgenome -T 8 >gemmap_{4}.log &
#trim nextera adapters
echo "Trimming adapters"
trim_galore --paired --nextera $2 $3

#Map with GEM mapper
echo "Aligning"
gem-mapper -p -q 'offset-33' -T 30 -I indexedgenome.gem -1 "${3}_1_val_1.fq" -2 "${3}_2_val_2.fq" -o "${3}/mapped_${3}_wt_GEM"
#output mapped_name_wt_GEM.map
echo "Converting gem to sam"
#gem-2-sam
gem-2-sam -I indexedgenome.gem -c -T 30 -i "${3}/mapped_${3}_wt_GEM.map" -o "${3}/mapped_${3}_wt_GEM.sam" -q 'offset-33' -l
echo "Shrinking sam file.."
shrinksam -v -i ${3}/"mapped_${3}_wt_GEM.sam" -k "${3}/mapped_${3}_wt_GEM_shrinked.sam"
echo "Correcting sam headers.."
python /mnt/bioinfonew/srikanth/Test_align/genomeindices/correct_sam.py "${3}/mapped_${3}_wt_GEM_shrinked.sam"
echo "Running pathoid"
pathoscope ID -alignFile "${3}/mapped_${3}_wt_GEM_shrinked_updated.sam" -fileType sam -expTag "${3}_GEM" -outDir ${3}
echo "File stored as ${3}_GEM-sam-report.tsv"

