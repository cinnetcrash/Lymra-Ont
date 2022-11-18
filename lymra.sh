!#/bin/bash


# This is the main script file, do not change it if you don't know what you are doing.
# Gültekin Ünal, 2022


# Enter kraken database location
# Demunizer

# QC Step
mkdir QC_RAW_READS
fastqc data/* -o QC_RAW_READS

myreads=*.fastq
RAWREADS_2=*2.fastq
mv $myreads 01_raw
cd 01_raw
for i in *fastq; do echo "${i}: $(grep "@" $i | wc -l) reads"; done


# Read Trimming for Nanopore
mkdir trimmed_reads
porechop -i ${1}.fastq -o trimmed_reads/$1.fastq --format fastq -t 4


# FASTQ Stats

python3 


#spades.py --careful -o $1_SPADES_OUT -1 trimmed_reads/$1.fastq 
#CHECK THIS CODE!


# Mapping all reads to the reference genome
minimap2 -x map-ont -d coronovirus.mmi reference_genom.fasta
samtools view -bS sample1.sam > sample1.bam
samtools sort sample1.bam -o sample1.sorted.bam
samtools index sample1.sorted.bam
samtools view -q 10 -c sample1.sorted.bam
bcftools mpileup -f sars-cov-2_genome.fasta sample1.sorted.marked.bam | bcftools call -mv -Ob -o sample1.bcf
freebayes –f sars-cov-2_genome.fasta sample1.sorted.marked.bam > sample1.vcf
vcftools --vvcf sample1.vcf --minQ 20 --recode --recode-INFO-all --out sample1_q20
bcftools view -i '%QUAL>=2' sample1.vcf
lofreq call -f sars-cov-2_genome.fasta -o sample1.vcf sample1.sorted.marked.bam 
samtools mpileup -uf reference.fasta mapped_reads.bam | bcftools view -cg - | vcfutils vcf2fq  > consensus.fasta
samtools mpileup -aa -A -d 0 -Q 0 <input.bam> | ivar consensus -p <prefix> 

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

busco -f -c 20 -m genome -i Desktop/Monkeypox_Paper/Polished/002_BeratK_homopolished.fasta -o 002_busco --auto-lineage-prok

