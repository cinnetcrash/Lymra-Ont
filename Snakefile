import os
import shutil
import pandas as pd
from glob import glob


configfile: 'config.yaml'
os.makedirs("data/samples", exist_ok=True)

samples = pd.read_csv('samples.csv', sep='\t')


def get_samples():
    return list(samples['sample'].unique())


def get_fastq(wildcards):
    fastqs = samples.query("sample=='{}'".format(wildcards.sample))[["fastq_folder"]].iloc[0]
    ret_str = f"{fastqs.fastq_folder}"
    return fastqs

rule all:
    input:
        expand("data/return_files/{sample}_report.html",sample=get_samples())

rule copy_fastqs_to_path:
    input:
        get_fastq
    output:
        sample="data/fastq/{sample}.fastq"
    conda:
        srcdir("envs/ncov.yml")
    shell:
        """
        artic guppyplex --directory {input} --output {output.sample}
        """

rule porechop_trim:
    input:
        read_file="data/samples/{sample}.fastq.gz"
    output:
        trimmed_output="trimmed/{sample}.trimmed.fastq"
    conda: 
        "envs/conda-porechop.yaml"
    shell: 
        "porechop -i {input.read_file} -o {output.trimmed_output} -t 8 --barcode_threshold 60 --barcode_diff 1"

rule kraken2_viral:
    input:
        trimmed_reads="trimmed/{sample}.trimmed.fastq",
        database="/../../media/cinnet/PortableSSD1/kraken2_db/human_genome/"
    output:
        report_kraken2="kraken2_reports/{sample}_viral.txt",
        classified_out="trimmed/{sample}_viral_reads.fastq"
    conda:
        "envs/conda-kraken2.yaml"
    shell:
        "kraken2 --db {input.database} {input.trimmed_reads} --report {output.report_kraken2} --classified-out {output.classified_out}"

rule flye:
    input:
        "trimmed/{sample}_viral_reads.fastq"
    output:
        "assembly/flye_{sample}/"
    conda:
        "envs/conda-flye.yaml"
    shell:
        "flye --nano-corr {input} --out-dir {output} --genome-size 0.2m --meta -t 8"

rule minimap2_sort:
    input:
        fastq="data/fastq/{sample}.fastq",
        reference=config["reference"]
    params:
        preset="-x map-ont", #Oxford Nanopore to reference mapping (-k15)
        threads=config["threads"],
        kmer="-k 21", #k-mer length
        secondary_aligment="--secondary=yes -N 5", #Whether to output secondary alignments with most INT secondary alignments
        sec_score = "-p 0.8" #Minimal secondary-to-primary score ratio to output secondary mappings
    output:
        bam="data/mapped/{sample}.sorted.bam",
        bai="data/mapped/{sample}.sorted.bam.bai"
    conda:
        "envs/flye.yaml"
    shell:
        """
        minimap2 -a --frag=yes -t {params.threads} {params.preset} {params.kmer} {params.secondary_aligment} \
        {params.sec_score} {input.reference} {input.fastq} | samtools view -bS - | samtools sort -o {output.bam} - ; \n
        samtools index {output.bam}
        """

rule medaka:
    input:
        fq="trimmed/{sample}_viral_reads.fastq",
        reference="data/monkeypox.fa"
    output:
        "medaka_output_{sample}/"
    conda:
        "envs/conda-medaka.yaml"
    shell:
        """
        medaka_consensus -i {input.fq} -d {input.reference} -o {output} -t 6 -m r941_min_high_g303
        mv medaka_output_{sample}/consensus.fasta medaka_output_{sample}/{sample}.fasta
        """

rule homopolish:
    conda: "envs/conda-homopolish.yaml"
    input:
        medaka_output_fasta="medaka_output_{sample}/{sample}.fasta",
        ref="data/monkeypox.fa"
    output: 
        fa="polish/homopolish/{sample}_polished.fasta"
    shell:
        "homopolish polish -a {input.medaka_output_fasta} -l {input.ref} -m R9.4.pkl -o {output.fa}"

rule busco:
    input:
        "polish/homopolish/{sample}_polished.fasta"
    output:
        "busco_output/{sample}_run/"
    conda:
        "envs/conda-porechop.yaml"
    shell:
        "busco -i {input} -o {output} --auto-lineage-prok -m genome"