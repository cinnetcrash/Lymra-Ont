!#/bin/bash

# This is the main script file, do not change it if you don't know what you are doing.
# Gültekin Ünal, 2022


# QC Step
mkdir QC_RAW_READS
fastqc data/* -o QC_RAW_READS

myreads=*.fastq
RAWREADS_2=*2.fastq
mv $RAWREADS_1 $RAWREADS_2 01_raw
cd 01_raw
for i in *fastq; do echo "${i}: $(grep "@" $i | wc -l) reads"; done


# Read Trimming for Nanopore
mkdir trimmed_reads
porechop -i ${1}.fastq -o trimmed_reads/$1.fastq --format fastq -t 4

#spades.py --careful -o $1_SPADES_OUT -1 trimmed_reads/$1.fastq 
#CHECK THIS CODE!


# Mapping all reads to the reference genome
minimap2


# Polishing draft assembly with Homopolish 

mkdir polished_assembly

#increase the number of threads if you have more cpu cores to make it run faster
threads=4

currdir=$(pwd)
pilon=$currdir/apps/pilon.jar
read1=$currdir/trimmed_reads/P7741_R1.fastq.gz
read2=$currdir/trimmed_reads/P7741_R2.fastq.gz  


mkdir polishing_process
cd polishing_process

echo "POLISING DRAFT ASSEMBLY"

# Homopolish Step

echo "All processes done"


