#!/bin/bash
#PBS -l select=5:ncpus=24:mpiprocs=24
#PBS -P CBBI1661
#PBS -q normal
#PBS -l walltime=48:00:00
#PBS -o /mnt/lustre/users/rhenriques/orange_roughy/pipeline1.out
#PBS -e /mnt/lustre/users/rhenriques/orange_roughy/pipeline2.err
#PBS -m abe
ulimit -s unlimited


cd /mnt/lustre/users/rhenriques/orange_roughy/03_bwa
nproc=`cat $PBS_NODEFILE | wc -l`
echo nproc is $nproc
cat $PBS_NODEFILE


module add chpc/BIOMODULES
module load fastqc/0.11.9
module load cutadapt/3.4
module load trim_galore/0.6.5
module load trimmomatic/0.36
module load bwa/0.7.17
module load samtools/1.9
module load stacks/2.66
module load bcftools/1.9
module load vcftools/0.1.15


# prepare the samples, as STACKS does not recognize .bz2 extension

 # Directory containing the bz2 files 
 directory="/mnt/lustre/users/rhenriques/orange_roughy/RAW/bzp2_folder" 
 # Loop through each bz2 file in the directory and unzip it 
 for file in "$directory"/*.bz2; do [ -f "$file" ] && echo "Unzipping $file" && bzip2 -d "$file" || echo "No bz2 files found in the directory."; done; echo "All files unzipped."



#run fastqc and multiqc
fastqc orange_roughy/RAW/*fastq -o orange_roughy/00_fastqc_multiqc
multiqc orange_roughy/00_fastqc_multiqc -o orange_roughy/00_fastqc_multiqc


# process radtags
#this one is working well now - no absolute paths, navigate to where the samples are and run process_radtags from there

cd orange_roughy/RAW

for i in $(cat /mnt/lustre/users/rhenriques/orange_roughy/samples_list2.txt); do process_radtags -P -i fastq -1 $i"_R1.fastq" -2 $i"_R2.fastq"  -o ../01_process_radtags --renz-1 pstI --renz-2 mspI -c -q -r --threads 4; done



# clean and trim reads to the same length

cd ..
for i in $(cat /mnt/lustre/users/rhenriques/orange_roughy/samples_list2.txt); do trim_galore -q 28 --length 50 --clip_R1 5 --clip_R2 5 --paired --gzip -o 02_trimmed/ 01_process_radtags/$i"_R1.1.fq" 01_process_radtags/$i"_R2.2.fq"; done

for i in $(cat /mnt/lustre/users/rhenriques/orange_roughy/samples_list2.txt); do trimmomatic PE -phred33 02_trimmed/$i"_R1.1_val_1.fq.gz" 02_trimmed/$i"_R2.2_val_2.fq.gz" 02_trimmed/$i"_1_paired.fq.gz" 02_trimmed/$i"_1_unpaired.fq.gz" 02_trimmed/$i"_2_paired.fq.gz" 02_trimmed/$i"_2_unpaired.fq.gz" MINLEN:135 CROP:135;done


# map  

# unzip the reference first, as samtools doesn't work with zipped files

gunzip orange_roughy_genome/*fna.gz

bwa index orange_roughy_genome/*fna
samtools faidx orange_roughy_genome/*fna

/apps/chpc/bio/picard/3.0.0/picard.0 CreateSequenceDictionary REFERENCE=orange_roughy_genome/*fna OUTPUT=orange_roughy_genome/orange_roughy.dict

for i in $(cat /mnt/lustre/users/rhenriques/orange_roughy/samples_list2.txt); do bwa mem -t 24 -M orange_roughy_genome/*fna 02_trimmed/$i"_1_paired.fq.gz" 02_trimmed/$i"_2_paired.fq.gz" > 03_bwa/$i".sam"; samtools view -@ 24 -Sb -q 20 -o 03_bwa/$i".bam" 03_bwa/$i".sam"; done

rm 03_bwa/*.sam 

for i in $(cat /mnt/lustre/users/rhenriques/orange_roughy/samples_list2.txt); do samtools sort 03_bwa/$i".bam" > 03_bwa/$i"_sorted.bam"; samtools index 03_bwa/$i"_sorted.bam"; done

for i in $(cat /mnt/lustre/users/rhenriques/orange_roughy/samples_list2.txt); do /apps/chpc/bio/picard/3.0.0/picard.0 AddOrReplaceReadGroups I=03_bwa/$i"_sorted.bam" O=03_bwa/$i"_RG.bam" RGID=$i RGLB=libOR RGPL=IlluminaNextSeq500 RGPU=Oct2023 RGSM=$i; samtools sort 03_bwa/$i"_RG.bam" > 03_bwa/$i"_sorted_RG.bam"; samtools index 03_bwa/$i"_sorted_RG.bam"; done


#keep only properly mapped and paired reads
for i in $(cat /mnt/lustre/users/rhenriques/orange_roughy/samples_list2.txt); do samtools view -b -f 2 -o 03_bwa/$i"_sorted_RG_primary.bam" 03_bwa/$i"_sorted_RG.bam"; done


#check if there are still secondary reads
for i in $(cat /mnt/lustre/users/rhenriques/orange_roughy/samples_list2.txt); do samtools flagstat 03_bwa/$i"_sorted_RG_primary2.bam" > 03_bwa/$i"_stats_primary2.txt"; done

for i in $(cat /mnt/lustre/users/rhenriques/orange_roughy/samples_list2.txt); do samtools view -b -F 256 -o 03_bwa/$i"_sorted_RG_primary2.bam" 03_bwa/$i"_sorted_RG_primary.bam"; done

for i in $(cat /mnt/lustre/users/rhenriques/orange_roughy/samples_list2.txt); do samtools sort 03_bwa/$i"_sorted_RG_primary2.bam" > 03_bwa/$i"_sorted_RG_primary2_sorted.bam"; samtools index 03_bwa/$i"_sorted_RG_primary2_sorted.bam"; done


for i in $(cat /mnt/lustre/users/rhenriques/orange_roughy/samples_list2.txt); do rm 03_bwa/$i".bam"; done


# call SNPs

# bcftools

cd 03_bwa

bcftools  mpileup --threads 120 -E -C50 -f ../orange_roughy_genome/*fna P01-A10-Plate1-R04_sorted_RG.bam P01-A11-Plate1-R14_sorted_RG.bam P01-B10-Plate1-R05_sorted_RG.bam P01-B11-Plate1-R15_sorted_RG.bam P01-C10-Plate1-R06_sorted_RG.bam P01-C11-Plate1-R16_sorted_RG.bam P01-D10-Plate1-R08_sorted_RG.bam P01-D11-Plate1-R17_sorted_RG.bam P01-E10-Plate1-R09_sorted_RG.bam P01-E11-Plate1-R18_sorted_RG.bam P01-F09-Plate1-R01_sorted_RG.bam P01-F10-Plate1-R10_sorted_RG.bam P01-F11-Plate1-R19_sorted_RG.bam P01-G09-Plate1-R02_sorted_RG.bam P01-G10-Plate1-R12_sorted_RG.bam P01-G11-Plate1-R20_sorted_RG.bam P01-H09-Plate1-R03_sorted_RG.bam P01-H10-Plate1-R13_sorted_RG.bam  P01-A08-Plate1-T08_sorted_RG.bam P01-B07-Plate1-T01_sorted_RG.bam P01-B08-Plate1-T09_sorted_RG.bam P01-C07-Plate1-T02_sorted_RG.bam P01-C08-Plate1-T10_sorted_RG.bam P01-D07-Plate1-T03_sorted_RG.bam P01-E07-Plate1-T04_sorted_RG.bam P01-F07-Plate1-T05_sorted_RG.bam P01-G07-Plate1-T06_sorted_RG.bam P01-H07-Plate1-T07_sorted_RG.bam P01-A09-Plate1-F06_sorted_RG.bam P01-B09-Plate1-F07_sorted_RG.bam P01-C09-Plate1-F08_sorted_RG.bam P01-D08-Plate1-F01_sorted_RG.bam P01-D09-Plate1-F09_sorted_RG.bam P01-E08-Plate1-F02_sorted_RG.bam P01-E09-Plate1-F10_sorted_RG.bam P01-F08-Plate1-F03_sorted_RG.bam P01-G08-Plate1-F04_sorted_RG.bam P01-H08-Plate1-F05_sorted_RG.bam  P02-A01-Plate2-H01_sorted_RG.bam P02-A02-Plate2-H09_sorted_RG.bam P02-A03-Plate2-H17_sorted_RG.bam P02-B01-Plate2-H02_sorted_RG.bam P02-B02-Plate2-H10_sorted_RG.bam P02-B03-Plate2-H18_sorted_RG.bam P02-C01-Plate2-H03_sorted_RG.bam P02-C02-Plate2-H11_sorted_RG.bam P02-C03-Plate2-H19_sorted_RG.bam P02-D01-Plate2-H04_sorted_RG.bam P02-D02-Plate2-H12_sorted_RG.bam P02-D03-Plate2-H20_sorted_RG.bam P02-E01-Plate2-H05_sorted_RG.bam P02-E02-Plate2-H13_sorted_RG.bam P02-F01-Plate2-H06_sorted_RG.bam P02-F02-Plate2-H14_sorted_RG.bam P02-G01-Plate2-H07_sorted_RG.bam P02-G02-Plate2-H15_sorted_RG.bam P02-H01-Plate2-H08_sorted_RG.bam P02-H02-Plate2-H16_sorted_RG.bam P02-A04-Plate2-J05_sorted_RG.bam P02-A05-Plate2-J13_sorted_RG.bam P02-B04-Plate2-J06_sorted_RG.bam P02-B05-Plate2-J14_sorted_RG.bam P02-C04-Plate2-J07_sorted_RG.bam P02-C05-Plate2-J15_sorted_RG.bam P02-D04-Plate2-J08_sorted_RG.bam P02-D05-Plate2-J16_sorted_RG.bam P02-E03-Plate2-J01_sorted_RG.bam P02-E04-Plate2-J09_sorted_RG.bam P02-E05-Plate2-J17_sorted_RG.bam P02-F03-Plate2-J02_sorted_RG.bam P02-F04-Plate2-J10_sorted_RG.bam P02-F05-Plate2-J18_sorted_RG.bam P02-G03-Plate2-J03_sorted_RG.bam P02-G04-Plate2-J11_sorted_RG.bam P02-G05-Plate2-J19_sorted_RG.bam P02-H03-Plate2-J04_sorted_RG.bam P02-H04-Plate2-J12_sorted_RG.bam P02-H05-Plate2-J20_sorted_RG.bam  P01-D03-Plate1-S01_sorted_RG.bam P01-E03-Plate1-S03_sorted_RG.bam P01-F03-Plate1-S04_sorted_RG.bam P01-G03-Plate1-S05_sorted_RG.bam P01-H03-Plate1-S09_sorted_RG.bam P01-A04-Plate1-S11_sorted_RG.bam P01-B04-Plate1-S28_sorted_RG.bam P01-C04-Plate1-S29_sorted_RG.bam P01-D04-Plate1-S30_sorted_RG.bam P01-E04-Plate1-S31_sorted_RG.bam P01-F04-Plate1-S32_sorted_RG.bam P01-G04-Plate1-S35_sorted_RG.bam P01-H04-Plate1-S39_sorted_RG.bam P01-A05-Plate1-S40_sorted_RG.bam P01-B05-Plate1-S41_sorted_RG.bam P01-C05-Plate1-S42_sorted_RG.bam P01-D05-Plate1-S44_sorted_RG.bam P01-E05-Plate1-S45_sorted_RG.bam P01-F05-Plate1-S46_sorted_RG.bam P01-G05-Plate1-S47_sorted_RG.bam P01-H05-Plate1-S48_sorted_RG.bam P01-A06-Plate1-S49_sorted_RG.bam P01-B06-Plate1-S51_sorted_RG.bam P01-C06-Plate1-S52_sorted_RG.bam P01-D06-Plate1-S53_sorted_RG.bam P01-E06-Plate1-S54_sorted_RG.bam P01-F06-Plate1-S55_sorted_RG.bam P01-G06-Plate1-S56_sorted_RG.bam P01-H06-Plate1-S57_sorted_RG.bam P01-A07-Plate1-S58_sorted_RG.bam  P01-A01-Plate1-Z01_sorted_RG.bam P01-A02-Plate1-Z10_sorted_RG.bam P01-A03-Plate1-Z21_sorted_RG.bam P01-B01-Plate1-Z02_sorted_RG.bam P01-B02-Plate1-Z11_sorted_RG.bam P01-B03-Plate1-Z22_sorted_RG.bam P01-C01-Plate1-Z03_sorted_RG.bam P01-C02-Plate1-Z12_sorted_RG.bam P01-C03-Plate1-Z23_sorted_RG.bam P01-D01-Plate1-Z04_sorted_RG.bam P01-D02-Plate1-Z13_sorted_RG.bam P01-E01-Plate1-Z05_sorted_RG.bam P01-E02-Plate1-Z15_sorted_RG.bam P01-F01-Plate1-Z07_sorted_RG.bam P01-F02-Plate1-Z17_sorted_RG.bam P01-G01-Plate1-Z08_sorted_RG.bam P01-G02-Plate1-Z19_sorted_RG.bam P01-H01-Plate1-Z09_sorted_RG.bam P01-H02-Plate1-Z20_sorted_RG.bam | bcftools call --threads 120 -mv -Oz -o ../04_snp_calling_bcftools/all_bcftools3.vcf.gz

bcftools reheader -s  ../04_snp_calling_bcftools/all_bcftools_filtered.vcf -o ../04_snp_calling_bcftools/all_bcftools_filtered.vcf

#bcftools view --threads 24 -O z -e 'F_MISSING>0.0 || MAF<0.05 || DP<5 || DP > 100' --types snps -m2 -M2 ../04_snp_calling_bcftools/all_bcftools.vcf -o ../04_snp_calling_bcftools/all_bcftools_filtered.vcf
 
/apps/chpc/bio/freebayes/1.3.1/bin/freebayes -f ../orange_roughy_genome/*fna --min-coverage 3 -F 0.05 P01-A10-Plate1-R04_sorted_RG.bam P01-A11-Plate1-R14_sorted_RG.bam P01-B10-Plate1-R05_sorted_RG.bam P01-B11-Plate1-R15_sorted_RG.bam P01-C10-Plate1-R06_sorted_RG.bam P01-C11-Plate1-R16_sorted_RG.bam P01-D10-Plate1-R08_sorted_RG.bam P01-D11-Plate1-R17_sorted_RG.bam P01-E10-Plate1-R09_sorted_RG.bam P01-E11-Plate1-R18_sorted_RG.bam P01-F09-Plate1-R01_sorted_RG.bam P01-F10-Plate1-R10_sorted_RG.bam P01-F11-Plate1-R19_sorted_RG.bam P01-G09-Plate1-R02_sorted_RG.bam P01-G10-Plate1-R12_sorted_RG.bam P01-G11-Plate1-R20_sorted_RG.bam P01-H09-Plate1-R03_sorted_RG.bam P01-H10-Plate1-R13_sorted_RG.bam  P01-A08-Plate1-T08_sorted_RG.bam P01-B07-Plate1-T01_sorted_RG.bam P01-B08-Plate1-T09_sorted_RG.bam P01-C07-Plate1-T02_sorted_RG.bam P01-C08-Plate1-T10_sorted_RG.bam P01-D07-Plate1-T03_sorted_RG.bam P01-E07-Plate1-T04_sorted_RG.bam P01-F07-Plate1-T05_sorted_RG.bam P01-G07-Plate1-T06_sorted_RG.bam P01-H07-Plate1-T07_sorted_RG.bam P01-A09-Plate1-F06_sorted_RG.bam P01-B09-Plate1-F07_sorted_RG.bam P01-C09-Plate1-F08_sorted_RG.bam P01-D08-Plate1-F01_sorted_RG.bam P01-D09-Plate1-F09_sorted_RG.bam P01-E08-Plate1-F02_sorted_RG.bam P01-E09-Plate1-F10_sorted_RG.bam P01-F08-Plate1-F03_sorted_RG.bam P01-G08-Plate1-F04_sorted_RG.bam P01-H08-Plate1-F05_sorted_RG.bam  P02-A01-Plate2-H01_sorted_RG.bam P02-A02-Plate2-H09_sorted_RG.bam P02-A03-Plate2-H17_sorted_RG.bam P02-B01-Plate2-H02_sorted_RG.bam P02-B02-Plate2-H10_sorted_RG.bam P02-B03-Plate2-H18_sorted_RG.bam P02-C01-Plate2-H03_sorted_RG.bam P02-C02-Plate2-H11_sorted_RG.bam P02-C03-Plate2-H19_sorted_RG.bam P02-D01-Plate2-H04_sorted_RG.bam P02-D02-Plate2-H12_sorted_RG.bam P02-D03-Plate2-H20_sorted_RG.bam P02-E01-Plate2-H05_sorted_RG.bam P02-E02-Plate2-H13_sorted_RG.bam P02-F01-Plate2-H06_sorted_RG.bam P02-F02-Plate2-H14_sorted_RG.bam P02-G01-Plate2-H07_sorted_RG.bam P02-G02-Plate2-H15_sorted_RG.bam P02-H01-Plate2-H08_sorted_RG.bam P02-H02-Plate2-H16_sorted_RG.bam P02-A04-Plate2-J05_sorted_RG.bam P02-A05-Plate2-J13_sorted_RG.bam P02-B04-Plate2-J06_sorted_RG.bam P02-B05-Plate2-J14_sorted_RG.bam P02-C04-Plate2-J07_sorted_RG.bam P02-C05-Plate2-J15_sorted_RG.bam P02-D04-Plate2-J08_sorted_RG.bam P02-D05-Plate2-J16_sorted_RG.bam P02-E03-Plate2-J01_sorted_RG.bam P02-E04-Plate2-J09_sorted_RG.bam P02-E05-Plate2-J17_sorted_RG.bam P02-F03-Plate2-J02_sorted_RG.bam P02-F04-Plate2-J10_sorted_RG.bam P02-F05-Plate2-J18_sorted_RG.bam P02-G03-Plate2-J03_sorted_RG.bam P02-G04-Plate2-J11_sorted_RG.bam P02-G05-Plate2-J19_sorted_RG.bam P02-H03-Plate2-J04_sorted_RG.bam P02-H04-Plate2-J12_sorted_RG.bam P02-H05-Plate2-J20_sorted_RG.bam  P01-D03-Plate1-S01_sorted_RG.bam P01-E03-Plate1-S03_sorted_RG.bam P01-F03-Plate1-S04_sorted_RG.bam P01-G03-Plate1-S05_sorted_RG.bam P01-H03-Plate1-S09_sorted_RG.bam P01-A04-Plate1-S11_sorted_RG.bam P01-B04-Plate1-S28_sorted_RG.bam P01-C04-Plate1-S29_sorted_RG.bam P01-D04-Plate1-S30_sorted_RG.bam P01-E04-Plate1-S31_sorted_RG.bam P01-F04-Plate1-S32_sorted_RG.bam P01-G04-Plate1-S35_sorted_RG.bam P01-H04-Plate1-S39_sorted_RG.bam P01-A05-Plate1-S40_sorted_RG.bam P01-B05-Plate1-S41_sorted_RG.bam P01-C05-Plate1-S42_sorted_RG.bam P01-D05-Plate1-S44_sorted_RG.bam P01-E05-Plate1-S45_sorted_RG.bam P01-F05-Plate1-S46_sorted_RG.bam P01-G05-Plate1-S47_sorted_RG.bam P01-H05-Plate1-S48_sorted_RG.bam P01-A06-Plate1-S49_sorted_RG.bam P01-B06-Plate1-S51_sorted_RG.bam P01-C06-Plate1-S52_sorted_RG.bam P01-D06-Plate1-S53_sorted_RG.bam P01-E06-Plate1-S54_sorted_RG.bam P01-F06-Plate1-S55_sorted_RG.bam P01-G06-Plate1-S56_sorted_RG.bam P01-H06-Plate1-S57_sorted_RG.bam P01-A07-Plate1-S58_sorted_RG.bam  P01-A01-Plate1-Z01_sorted_RG.bam P01-A02-Plate1-Z10_sorted_RG.bam P01-A03-Plate1-Z21_sorted_RG.bam P01-B01-Plate1-Z02_sorted_RG.bam P01-B02-Plate1-Z11_sorted_RG.bam P01-B03-Plate1-Z22_sorted_RG.bam P01-C01-Plate1-Z03_sorted_RG.bam P01-C02-Plate1-Z12_sorted_RG.bam P01-C03-Plate1-Z23_sorted_RG.bam P01-D01-Plate1-Z04_sorted_RG.bam P01-D02-Plate1-Z13_sorted_RG.bam P01-E01-Plate1-Z05_sorted_RG.bam P01-E02-Plate1-Z15_sorted_RG.bam P01-F01-Plate1-Z07_sorted_RG.bam P01-F02-Plate1-Z17_sorted_RG.bam P01-G01-Plate1-Z08_sorted_RG.bam P01-G02-Plate1-Z19_sorted_RG.bam P01-H01-Plate1-Z09_sorted_RG.bam P01-H02-Plate1-Z20_sorted_RG.bam --vcf ../04_snp_calling_bcftools/all_freebayes.vcf

# stacks

denovo_map.pl -T 8 --samples $samples --popmap $pop_map1 --out-path $denovo --paired --min-samples-per-pop 0.80 -M 2 -m 3 -n 2 -X "populations: --min-maf 0.05 --write-random-snp --vcf --fstats --hwe"
denovo_map.pl -T 8 --samples $samples --popmap $pop_map1 --out-path $denovo --paired --min-samples-per-pop 0.80 -M 2 -m 3 -n 2 -X "populations: --min-maf 0.05 --write-random-snp --vcf --fstats --hwe"


# filter SNPs

#note, if it says primary2/prim, it means the bams were filtered for only properly paired and no secondary alignments
vcftools --gzvcf 04_snp_calling_bcftools/all_bcftools_primary2_2.vcf.gz --min-alleles 2 --max-alleles 2 --recode --recode-INFO-all --out 04_snp_calling_bcftools/1.bcftools.primary2_2
vcftools --vcf 04_snp_calling_bcftools/1.bcftools.primary2_2.recode.vcf --remove-indels --recode --recode-INFO-all --out 04_snp_calling_bcftools/2.bcftools.primary2_2

vcftools --vcf 04_snp_calling_bcftools/2.bcftools.primary2_2.recode.vcf --missing-site --out 04_snp_calling_bcftools/2.MissSite_prim2
cat 04_snp_calling_bcftools/2.MissSite_prim2.lmiss | awk '!/CHR/' | awk '$6 > 0.00' | cut -f1,2 >> 04_snp_calling_bcftools/2.badloci.prim2
  vcftools --vcf 04_snp_calling_bcftools/2.bcftools.primary2_2.recode.vcf --exclude-positions 04_snp_calling_bcftools/2.badloci.prim2 --recode --recode-INFO-all --out 04_snp_calling_bcftools/3.bcftools.primary2_2

vcftools --vcf 04_snp_calling_bcftools/5.bcftools.primary2_2.recode.vcf --missing-indv --out 04_snp_calling_bcftools/5.MissIndv.prim2
cat 04_snp_calling_bcftools/3.MissIndv.prim2.imiss | awk '!/CHR/' | awk '$5 > 0.05' | cut -f1,2 >> 04_snp_calling_bcftools/3.badindv.prim2
vcftools --vcf 04_snp_calling_bcftools/3.bcftools3.recode.vcf --remove 04_snp_calling_bcftools/3.badindv --recode --recode-INFO-all --out 04_snp_calling_bcftools/4.bcftools3

vcftools --vcf 04_snp_calling_bcftools/4.bcftools3.recode.vcf --missing-site --out 04_snp_calling_bcftools/4.MissSite
cat 04_snp_calling_bcftools/4.MissSite.lmiss | awk '!/CHR/' | awk '$6 > 0.05' | cut -f1,2 >> 04_snp_calling_bcftools/4.badloci

vcftools --vcf 04_snp_calling_bcftools/3.bcftools.primary2_2.recode.vcf --mac 3 --recode --recode-INFO-all --out 04_snp_calling_bcftools/4.bcftools.primary2_2

vcftools --vcf 04_snp_calling_bcftools/4.bcftools.primary2_2.recode.vcf --minDP 5 --recode --recode-INFO-all --out 04_snp_calling_bcftools/5.bcftools.primary2_2

vcftools --vcf 04_snp_calling_bcftools/5.bcftools.primary2.recode.vcf --missing-site --out 04_snp_calling_bcftools/5.MissSite2_prim
cat 04_snp_calling_bcftools/5.MissSite_prim2.lmiss | awk '!/CHR/' | awk '$6 > 0.05' | cut -f1,2 >> 04_snp_calling_bcftools/5.badloci.prim2_005
  vcftools --vcf 04_snp_calling_bcftools/5.bcftools.primary2_2.recode.vcf --exclude-positions 04_snp_calling_bcftools/5.badloci.prim2 --recode --recode-INFO-all --out 04_snp_calling_bcftools/6.bcftools.primary2_2

vcftools --vcf 04_snp_calling_bcftools/5.bcftools.primary2_2.recode.vcf --hardy --out 04_snp_calling_bcftools/5.hwe2

cat 04_snp_calling_bcftools/5.hwe2.hwe | cut -f1,2,3 | awk '{split($3, counts, "/"); print $1, $2, counts[1], counts[2], counts[3]}' | awk '{print $1, $2, $3, $4, $5, ($4 / 126)}' | awk '$6 >= 0.9' > 04_snp_calling_bcftools/5.Het_Sites_prim2

cat 04_snp_calling_bcftools/5.Het_Sites_prim2 | cut -f1,2 >> 04_snp_calling_bcftools/5.badloci2_prim2

vcftools --vcf 04_snp_calling_bcftools/5.bcftools.primary2_2.recode.vcf --exclude-positions 04_snp_calling_bcftools/5.badloci2_prim2 --recode --recode-INFO-all --out 04_snp_calling_bcftools/6.bcftools.primary2_2

vcftools --vcf 04_snp_calling_bcftools/6.bcftools.primary2_2.recode.vcf --missing-site --out 04_snp_calling_bcftools/6.MissSite2.prim2
cat 04_snp_calling_bcftools/6.MissSite.prim2.lmiss | awk '!/CHR/' | awk '$6 > 0.05' | cut -f1,2 >> 04_snp_calling_bcftools/6.badloci2.prim2


vcftools --vcf 04_snp_calling_bcftools/6.bcftools.primary2.recode.vcf --exclude-positions 04_snp_calling_bcftools/6.badloci2.prim --recode --recode-INFO-all --out 04_snp_calling_bcftools/7.bcftools.primary2
#note, 7.bcftools3.recode.vcf has loci with 5% missing data, 7.2.bcftools.recode.vcf has loci with 0% missing data


vcftools --vcf 04_snp_calling_bcftools/7.bcftools.primary2_2.recode.vcf --missing-indv --out 04_snp_calling_bcftools/7.MissIndv.prim2


vcftools --vcf 04_snp_calling_bcftools/7.bcftools.primary2_2.recode.vcf --maf 0.05 --recode --recode-INFO-all --out 04_snp_calling_bcftools/8.bcftools.primary2_2
bcftools +prune -l 0.2 -w 1000 04_snp_calling_bcftools/8.bcftools.primary2_2.recode.vcf -Ov -o 04_snp_calling_bcftools/9.bcftools.primary2_2.vcf

vcftools --vcf 04_snp_calling_bcftools/9.bcftools.primary2_2.vcf --relatedness2 --out 04_snp_calling_bcftools/9.relate2_2
vcftools --vcf 04_snp_calling_bcftools/9.bcftools.primary2.vcf --remove 04_snp_calling_bcftools/remove_related_primary2 --recode --recode-INFO-all --out 04_snp_calling_bcftools/10.bcftools.primary2


vcftools --vcf 04_snp_calling_bcftools/10.bcftools.primary2_2.recode.vcf --missing-site --out 04_snp_calling_bcftools/10.MissSite.prim2
cat 04_snp_calling_bcftools/10.MissSite.prim2.lmiss | awk '!/CHR/' | awk '$6 > 0.00' | cut -f1,2 >> 04_snp_calling_bcftools/10.badloci2.prim2
vcftools --vcf 04_snp_calling_bcftools/10.2.bcftools3.recode.vcf --exclude-positions 04_snp_calling_bcftools/10.2.badloci2 --recode --recode-INFO-all --out 04_snp_calling_bcftools/11.2.bcftools3
vcftools --vcf 04_snp_calling_bcftools/10.bcftools.primary2.recode.vcf --missing-indv --out 04_snp_calling_bcftools/10.MissIndv.prim


vcftools --vcf 04_snp_calling_bcftools/11.2.bcftools3.1.recode2.vcf --positions 04_snp_calling_bcftools/snps_pcadapt5 --recode --recode-INFO-all --out 04_snp_calling_bcftools/12.bcftools3

