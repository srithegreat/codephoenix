Adaptor trimming
----------------
java -jar /home/biouser1/tools/Trimmomatic-0.36/trimmomatic-0.36.jar PE -threads 52 -trimlog /dev/null -phred33 00_raw_data/SBRN18_R1.fastq.gz 00_raw_data/SBRN18_R2.fastq.gz 01_adapter_trimmed_data/AT_SBRN18_R1.fastq.gz /dev/null 01_adapter_trimmed_data/AT_SBRN18_R2.fastq.gz /dev/null ILLUMINACLIP:/home/biouser1/tools/Trimmomatic-0.36/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:15 MINLEN:20 1> 01_adapter_trimmed_data/AT_SBRN18.log 2> 01_adapter_trimmed_data/AT_SBRN18.err

java -jar /home/biouser1/tools/Trimmomatic-0.36/trimmomatic-0.36.jar PE -threads 52 -trimlog /dev/null -phred33 00_raw_data/SBRT18_R1.fastq.gz 00_raw_data/SBRT18_R2.fastq.gz 01_adapter_trimmed_data/AT_SBRT18_R1.fastq.gz /dev/null 01_adapter_trimmed_data/AT_SBRT18_R2.fastq.gz /dev/null ILLUMINACLIP:/home/biouser1/tools/Trimmomatic-0.36/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:15 MINLEN:20 1> 01_adapter_trimmed_data/AT_SBRT18.log 2> 01_adapter_trimmed_data/AT_SBRT18.err

FastQC
------
/home/biouser1/tools/FastQC/fastqc --nogroup -t 8 00_raw_data/*.gz 01_adapter_trimmed_data/*.gz

Alignment
---------
Template: nohup bwa mem -t 52 /home/data/bamReCalibration/HG19_BWA_INDEX/HG19_BWA_INDEX R1.fastq R2.fastq | samtools view -bS - | samtools sort -O BAM - > sample.sorted.bam 

bwa mem -R '@RG\tID:1\tSM:SBRN18\tPL:Illumina\tLB:lib1\tPU:unit1' -t 25 /home/data/bamReCalibration/HG19_BWA_INDEX/HG19_BWA_INDEX /home/data/exomeLocalPipeLine/01_adapter_trimmed_data/AT_SBRN18_R1.fastq.gz /home/data/exomeLocalPipeLine/01_adapter_trimmed_data/AT_SBRN18_R2.fastq.gz  | samtools view -bS - | samtools sort -O BAM -T temp -o SBRN18.sorted.bam -m 1G - |& tee SBRN18.sorted.err
bwa mem -R '@RG\tID:1\tSM:SBRT18\tPL:Illumina\tLB:lib1\tPU:unit1' -t 25 /home/data/bamReCalibration/HG19_BWA_INDEX/HG19_BWA_INDEX /home/data/exomeLocalPipeLine/01_adapter_trimmed_data/AT_SBRT18_R1.fastq.gz /home/data/exomeLocalPipeLine/01_adapter_trimmed_data/AT_SBRT18_R2.fastq.gz  | samtools view -bS - | samtools sort -O BAM -T temp -o SBRT18.sorted.bam -m 1G - |& tee SBRT18.sorted.err

samtools index sample.sorted.bam 

DeDup
-----
java -Xmx100G -Xms50G -jar /home/biouser1/tools/picard-tools-2.9.0/picard.jar MarkDuplicates REMOVE_DUPLICATES=FALSE INPUT=/home/biouser1/Desktop/data/exomeLocalPipeLine/02_sorted_bam/tumor/SBRT18.sorted.bam METRICS_FILE=SBRT18.metrics.txt OUTPUT=SBRT18.sorted.dedup.bam
java -Xmx100G -Xms50G -jar /home/biouser1/tools/picard-tools-2.9.0/picard.jar MarkDuplicates REMOVE_DUPLICATES=true INPUT=/home/biouser1/Desktop/data/exomeLocalPipeLine/02_sorted_bam/normal/SBRN18.sorted.bam METRICS_FILE=SBRN18.metrics.txt OUTPUT=SBRN18.sorted.dedup.bam

BQSR
----
java -Xmx100G -Xms50G -jar  $GATK -T BaseRecalibrator -R /home/data/bamReCalibration/HG19_BWA_INDEX/HG19_BWA_INDEX.fa -I /home/data/exomeLocalPipeLine/03_dedup_bam/SBRN18.sorted.dedup.bam -knownSites /home/data/bamReCalibration/knownIndels/1000G/1000G_phase1.indels.hg19.rehead.sorted.vcf -knownSites /home/data/bamReCalibration/knownIndels/Mills_and_1000G/Mills_and_1000G_gold_standard.indels.hg19.rehead.sorted.vcf -knownSites /home/data/bamReCalibration/dbSNP149/All_20161121.sorted.vcf -o SBRN18.sorted.dedup.grp -nct 52 -U ALLOW_SEQ_DICT_INCOMPATIBILITY 1> SBRN18.BaseRecalibrator.log 2> SBRN18.BaseRecalibrator.err &
java -Xmx100G -Xms50G -jar  $GATK -T BaseRecalibrator -R /home/data/bamReCalibration/HG19_BWA_INDEX/HG19_BWA_INDEX.fa -I /home/data/exomeLocalPipeLine/03_dedup_bam/SBRT18.sorted.dedup.bam -knownSites /home/data/bamReCalibration/knownIndels/1000G/1000G_phase1.indels.hg19.rehead.sorted.vcf -knownSites /home/data/bamReCalibration/knownIndels/Mills_and_1000G/Mills_and_1000G_gold_standard.indels.hg19.rehead.sorted.vcf -knownSites /home/data/bamReCalibration/dbSNP149/All_20161121.sorted.vcf -o SBRT18.sorted.dedup.grp -nct 52 -U ALLOW_SEQ_DICT_INCOMPATIBILITY 1> SBRT18.BaseRecalibrator.log 2> SBRT18.BaseRecalibrator.err &

IndelRealignment
----------------
java -Xmx100G -Xms50G -jar $GATK -T RealignerTargetCreator -R /home/data/bamReCalibration/HG19_BWA_INDEX/HG19_BWA_INDEX.fa -known /home/data/bamReCalibration/knownIndels/1000G/1000G_phase1.indels.hg19.rehead.sorted.vcf -known /home/data/bamReCalibration/knownIndels/Mills_and_1000G/Mills_and_1000G_gold_standard.indels.hg19.rehead.sorted.vcf -I ../04_BaseRecalibration/SBRN18.sorted.dedup.bqsr.bam -o SBRN18.sorted.dedup.bqsr.intervals -U ALLOW_SEQ_DICT_INCOMPATIBILITY 1> SBRN18.RealignerTargetCreator.log 2> SBRN18.RealignerTargetCreator.err &
java -Xmx100G -Xms50G -jar $GATK -T RealignerTargetCreator -R /home/data/bamReCalibration/HG19_BWA_INDEX/HG19_BWA_INDEX.fa -known /home/data/bamReCalibration/knownIndels/1000G/1000G_phase1.indels.hg19.rehead.sorted.vcf -known /home/data/bamReCalibration/knownIndels/Mills_and_1000G/Mills_and_1000G_gold_standard.indels.hg19.rehead.sorted.vcf -I ../04_BaseRecalibration/SBRT18.sorted.dedup.bqsr.bam -o SBRT18.sorted.dedup.bqsr.intervals -U ALLOW_SEQ_DICT_INCOMPATIBILITY 1> SBRT18.RealignerTargetCreator.log 2> SBRT18.RealignerTargetCreator.err &

java -Xmx100G -Xms50G -jar $GATK -T IndelRealigner -R /home/data/bamReCalibration/HG19_BWA_INDEX/HG19_BWA_INDEX.fa -I ../04_BaseRecalibration/SBRN18.sorted.dedup.bqsr.bam -o SBRN18.sorted.dedup.bqsr.indelRealigned.bam -targetIntervals SBRN18.sorted.dedup.bqsr.intervals -known /home/data/bamReCalibration/knownIndels/Mills_and_1000G/Mills_and_1000G_gold_standard.indels.hg19.rehead.sorted.vcf -known /home/data/bamReCalibration/knownIndels/1000G/1000G_phase1.indels.hg19.rehead.sorted.vcf 1> SBRN18.IndelRealigner.log 2> SBRN18.IndelRealigner.err &
java -Xmx100G -Xms50G -jar $GATK -T IndelRealigner -R /home/data/bamReCalibration/HG19_BWA_INDEX/HG19_BWA_INDEX.fa -I ../04_BaseRecalibration/SBRT18.sorted.dedup.bqsr.bam -o SBRT18.sorted.dedup.bqsr.indelRealigned.bam -targetIntervals SBRT18.sorted.dedup.bqsr.intervals -known /home/data/bamReCalibration/knownIndels/Mills_and_1000G/Mills_and_1000G_gold_standard.indels.hg19.rehead.sorted.vcf -known /home/data/bamReCalibration/knownIndels/1000G/1000G_phase1.indels.hg19.rehead.sorted.vcf 1> SBRT18.IndelRealigner.log 2> SBRT18.IndelRealigner.err &

MuTect2
-------
java -Xmx100G -Xms50G -jar $GATK -T MuTect2 -R /home/biouser1/Desktop/data/bamReCalibration/HG19_BWA_INDEX/HG19_BWA_INDEX.fa -I:normal /home/data/exomeLocalPipeLine/05_IndelRealignment/SBRN18.sorted.dedup.bqsr.indelRealigned.bam -I:tumor /home/data/exomeLocalPipeLine/05_IndelRealignment/SBRT18.sorted.dedup.bqsr.indelRealigned.bam --dbsnp /home/biouser1/Desktop/data/bamReCalibration/dbSNP149/All_20161121.sorted.vcf --cosmic /home/biouser1/Desktop/data/bamReCalibration/cosmic/CosmicCodingMuts.chrRenamed.sorted.vcf -o SBRT18.somatic.vcf -L /home/data/exomeLocalPipeLine/SureSelect_V5_50UTR/S04380219_Padded.bed -nct 25 -U ALLOW_SEQ_DICT_INCOMPATIBILITY 1> MuTect2.SBR_18.log 2> MuTect2.SBR_18.err &

MuTect1
-------
/home/biouser1/tools/jdk1.6.0_45/bin/java -Xmx100G -Xms50G -jar /home/biouser1/tools/MuTect1/muTect-1.1.4.jar --analysis_type MuTect --reference_sequence /home/biouser1/Desktop/data/bamReCalibration/HG19_BWA_INDEX/HG19_BWA_INDEX.fa --cosmic /home/biouser1/Desktop/data/bamReCalibration/cosmic/CosmicCodingMuts.chrRenamed.sorted.v4.1.vcf --dbsnp /home/biouser1/Desktop/data/bamReCalibration/dbSNP149/All_20161121.sorted.v4.1.vcf --intervals /home/data/exomeLocalPipeLine/SureSelect_V5_50UTR/S04380219_Padded.bed --input_file:normal /home/data/exomeLocalPipeLine/05_IndelRealignment/SBRN18.sorted.dedup.bqsr.indelRealigned.bam --input_file:tumor /home/data/exomeLocalPipeLine/05_IndelRealignment/SBRT18.sorted.dedup.bqsr.indelRealigned.bam --out SBRT18.somatic.vcf 1> MuTect1.SBR_18.log 2> MuTect1.SBR_18.err &

VCF downgrade
-------------
sed -e "s/##fileformat=VCFv4.2/##fileformat=VCFv4.1/"                All_20161121.sorted.vcf | sed "s/(//" | sed "s/)//" | sed "s/,Version=\"3\">/>/"  >                All_20161121.sorted.v4.1.vcf
sed -e "s/##fileformat=VCFv4.2/##fileformat=VCFv4.1/" CosmicCodingMuts.chrRenamed.sorted.vcf | sed "s/(//" | sed "s/)//" | sed "s/,Version=\"3\">/>/"  > CosmicCodingMuts.chrRenamed.sorted.v4.1.vcf

# If you are trying to view VCF 4.2 files in IGV - you may run into issues. This function might help you.
# This script will:
# 1. Rename the file as version 4.1
# 2. Replace parentheses in the INFO lines (IGV doesn't like these!)
# 
# function vcf_downgrade() {
#   outfile=${1/.bcf/}
#   outfile=${outfile/.gz/}
#   outfile=${outfile/.vcf/}
#   bcftools view --max-alleles 2 -O v $1 | \
#   sed "s/##fileformat=VCFv4.2/##fileformat=VCFv4.1/" | \
#   sed "s/(//" | \
#   sed "s/)//" | \
#   sed "s/,Version=\"3\">/>/" | \
#   bcftools view -O z > ${outfile}.dg.vcf.gz
#   tabix ${outfile}.dg.vcf.gz
# }

Annovar
-------
perl /home/biouser1/tools/annovar/table_annovar.pl SBRT18.somatic.vcf /home/biouser1/tools/annovar/humandb -buildver hg19 -out myanno -protocol refGene,snp138,dbnsfp33a,clinvar_20170130,icgc21,cosmic70,exac03 -operation g,f,f,f,f,f,f -nastring . -vcfinput

TODO
----
Indels from MuTect2 and SNPs from MuTect1 -> new VCF
Update recent version all BDs in Annovar; particularly dbSNP and Cosmic
Use convert2annovar.pl to convert VCF to Annovar format

