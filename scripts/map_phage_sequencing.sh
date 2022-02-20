#------------------------------------------------------------------#
# processing and variant calling of the re-sequenced phage data ####
#------------------------------------------------------------------#

# download reference genome and place it in ~/PRJEB50009/ref
# used https://www.ncbi.nlm.nih.gov/assembly/GCF_000886135.1
# Pseudomonas phage phi-2 (viruses)

mkdir -p ~/PRJEB50009/ref
# save as GCA_000886135.1_ViralProj42717_genomic.fna

#---------------------------#
# download data from ENA ####
#---------------------------#

# first download latest version of enaBrowserTools
wget https://github.com/enasequence/enaBrowserTools/archive/refs/tags/v0.0.3.zip
unzip v0.0.3.zip

# follow the instructions on enaBrowserTools GitHub
# https://github.com/enasequence/enaBrowserTools

# download all fastq files from project
enaGroupGet -f fastq -g read PRJEB50009

# change directory to project
cd PRJEB50009

ls

# ERR7976214 is phage 1, the resequenced lab strain of phi 2
# ERR7976218 is phage 2

#------------------#
# software used ####
#------------------#

# freebayes
# fastqc
# bcftools
# samtools
# minimap2
# samblaster
# vcflib
# trim-galore

#--------------------------#
# example code workflow ####
#--------------------------#

# first install mamba as its faster for installing than conda
conda install -c conda-forge mamba

# create an environment
mamba create -n sequencing_pipeline 

# open the environment
conda activate sequencing_pipeline

# install packages
conda install -c bioconda freebayes fastqc bcftools samtools minimap2 samblaster vcflib
conda install -c bioconda biopet-vcffilter
conda install -c bioconda trim-galore

#-------------------------------------#
# map and call variants of phage 1 ####
#-------------------------------------#

# list anaconda environments
conda env list

conda activate sequencing_pipeline

cd ERR7976214

# set working and other directories
wd=~/PRJEB50009/ERR7976214
trimmed_files=$wd/trimmed_files
fastqc_output=$wd/fastqc_output

mkdir -p $trimmed_files
mkdir -p $fastqc_output

# trim illumina reads
trim_galore --paired --quality 10 --output_dir $trimmed_files $wd/ERR7976214_1.fastq.gz $wd/ERR7976214_2.fastq.gz

# run fastqc on files
fastqc -o $fastqc_output $trimmed_files/*.fq.gz

# make directories for bam files
bam_files=$wd/bams
mkdir -p $bam_files

ref_phi2=~/PRJEB50009/ref/GCA_000886135.1_ViralProj42717_genomic.fna

# map reads to reference genome using minimap2
# use samblaster to mark duplicates
minimap2 -ax sr $ref_phi2 $trimmed_files/ERR7976214_1_val_1.fq.gz $trimmed_files/ERR7976214_2_val_2.fq.gz | samblaster | samtools view -S -b -o $bam_files/phage_1_minimap2.bam

# make directories for mapped and unmapped files
mkdir -p $bam_files/mapped
mkdir -p $bam_files/unmapped
mkdir -p $bam_files/stats

# 1. split bam into mapped and unmapped reads
samtools view -b -F4 $bam_files/phage_1_minimap2.bam > $bam_files/mapped/phage_1_mapped.bam
samtools view -f4 $bam_files/phage_1_minimap2.bam > $bam_files/unmapped/phage_1_unmapped.bam

# 2. sort mapped file by position in genome and not by order of mapped sequence
samtools sort -o $bam_files/mapped/phage_1_sorted.bam $bam_files/mapped/phage_1_mapped.bam

# 3. index the sorted bam file
samtools index $bam_files/mapped/phage_1_sorted.bam

# 4. remove intermediate (unsorted mapped bam file)
rm $bam_files/mapped/phage_1_mapped.bam

# 5. extract stats of file
samtools flagstat $bam_files/phage_1_minimap2.bam > $bam_files/stats/phage_1_raw_stats.txt

# run freebayes to call SNPs and indels
# make directory for output
mkdir -p $wd/vcf_output

# call SNPs using FreeBayes, set ploidy to 1
freebayes -f $ref_phi2  -p 1 $bam_files/mapped/phage_1_sorted.bam > $wd/vcf_output/phage_1.vcf
# keep only snps with a QUAL score > 20, based on recommendation of Erik Garrison
#https://www.biostars.org/p/71299/
vcffilter -f "QUAL > 20" $wd/vcf_output/phage_1.vcf > $wd/vcf_output/phage_1_filter.vcf

#-------------------------------------#
# map and call variants of phage 2 ####
#-------------------------------------#

# set working directory and create bam folder
wd=~/PRJEB50009/ERR7976218

mkdir -p $wd/bam_files

cd $wd

# use minimap2 to map long reads
minimap2 -ax map-pb $ref_phi2 ERR7976218.fastq.gz > $wd/bam_files/phage_2_minimap2.bam

# sort mapped file by position in genome and not by order of mapped sequence
samtools view -b -F 4 $wd/bam_files/phage_2_minimap2.bam > $wd/bam_files/phage_2_mapped.bam
samtools sort -o $wd/bam_files/phage_2_sorted.bam $wd/bam_files/phage_2_mapped.bam

# index the sorted bam file
samtools index $wd/bam_files/phage_2_sorted.bam

# get flagstats
samtools flagstat $wd/bam_files/phage_2_minimap2.bam > $wd/bam_files/phage_2_stats.txt

# make directory for output
mkdir -p $wd/vcf_output

# call SNPs using FreeBayes
freebayes -f $ref_phi2 -p 1 $wd/bam_files/phage_2_sorted.bam > $wd/vcf_output/phage_2.vcf
        
# keep only snps with a QUAL score > 20, based on recommendation of Erik Garrison #https://www.biostars.org/p/71299/
vcffilter -f "QUAL > 20" $wd/vcf_output/phage_2.vcf > $wd/vcf_output/phage_2_filter.vcf

# phage_1_filter.vcf and phage_2_filter.vcf become phage_1.vcf and phage_2.vcf in the GitHub repository