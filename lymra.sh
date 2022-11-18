!#/bin/bash


mkdir QC_RAW_READS
fastqc data/* -o QC_RAW_READS


mkdir trimmed_reads
porechop -i ${1}.fastq -o trimmed_reads/$1.fastq --format fastq -t 4


#spades.py --careful -o $1_SPADES_OUT -1 trimmed_reads/$1.fastq 
#CHECK THIS CODE!

./polish.sh
#polishing illumina draft assembly using pilon

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

echo "ROUND 1"

cp $currdir/P7741_SPADES_OUT/scaffolds.fasta ./raw_assembly.fasta
bwa index raw_assembly.fasta
bwa mem -t $threads raw_assembly.fasta $read1 $read2 | samtools view - -Sb |samtools sort - -@"threads" -o mapping1.sorted.bam
samtools index mapping1.sorted.bam
java -Xmx40G -jar  $pilon --genome raw_assembly.fasta --fix all --changes --frags mapping1.sorted.bam --threads $threads --output pilon_stage1|tee stage1.pilon


echo "Round2"

bwa index pilon_stage1.fasta
bwa mem -t $threads pilon_stage1.fasta $read1 $read2| samtools view - -Sb | samtools sort - -@"$threads" -o mapping2.sorted.bam
samtools index mapping2.sorted.bam
java -Xmx40G -jar $pilon --genome pilon_stage1.fasta --fix all --changes --frags mapping2.sorted.bam --threads $threads --output pilon_stage2 | tee stage2.pilon


echo "Round3"
bwa index pilon_stage2.fasta
bwa mem -t $threads pilon_stage2.fasta $read1 $read2| samtools view - -Sb | samtools sort - -@"$threads" -o mapping3.sorted.bam
samtools index mapping3.sorted.bam
java -Xmx40G -jar $pilon --genome pilon_stage2.fasta --fix all --changes --frags mapping3.sorted.bam --threads $threads --output pilon_stage3 | tee stage3.pilon



echo "Round4"
bwa index pilon_stage3.fasta
bwa mem -t $threads pilon_stage3.fasta $read1 $read2| samtools view - -Sb | samtools sort - -@"$threads" -o mapping4.sorted.bam
samtools index mapping4.sorted.bam
java -Xmx40G -jar $pilon --genome pilon_stage3.fasta --fix all --changes --frags mapping4.sorted.bam --threads $threads --output pilon_stage4 | tee stage4.pilon


cp pilon_stage4.fasta ../P7741.polished.fasta
cd ../
cat polishing_process/pilon_stage4.changes


./reorder_contigs.sh

#reference genome
ref=genomes/Liflandii.fasta

ragtag.py scaffold $ref P7741.polished.fasta -o P7741_reordered


#extract the reordered contig with a custom python script
#the scripts accept name of the ragtag file containing the reordered contigs and accession number for the reference genome
#accession number is found in the first line of the reference genome fasta file

python extract_reordered.py P7741_reordered/ragtag.scaffolds.fasta NC_020133.1


./annotate.sh

#Genome annotation using prokka

#if you have more cpus you can increase the number of cpus
cpus=4

prokka --cpus $cpus --kingdom Bacteria --locustag P7741 --outdir P7741_annotation --prefix P7741 --addgenes P7741.reordered.fasta
./get_pseudo.pl P7741_annotation/P7741.faa | tee P7741_annotation/P7741.pseudo.txt

./dendogram.sh

mkdir tmp_fastas
cp genomes/*.fasta tmp_fastas
cp P7741.reordered.fasta tmp_fastas/
mkdir dendogram
dRep compare dendogram -g tmp_fastas/*.fasta
rm -fr tmp_fastas

echo "dendogram generated"
echo "output: ./dendogram"


