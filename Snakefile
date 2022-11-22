from glob import glob
import pandas as pd
shell.executable("/bin/bash")

configfile: "config.yaml"


rule all:
    "busco_output/{sample}_busco.txt"




rule porechop_trim:
    input:
        "data/samples/{sample}.fastq"
    output:
        temp("trimmed/{sample}.trimmed.fastq")
    conda: 
        "envs/porechop.yaml"
    shell: 
        "porechop -i {input} -o {output} -t 3"

rule index_reference:
    input:
        "data/monkeypox.fa"
    output:
        temp("mapped/reference/monkeypox.mmi")
    conda:
        "envs/mapping.yaml"
    shell:
        "minimap2 -x map-ont -d {output} {input} | cp data/monkeypox.fa mapped/reference/"

rule flye:
    input:
        "data/{sample}.fastq"
    output:
        "assembly/{sample}.fasta"
    conda:
        "envs/flye.yaml"
    shell:
        "flye --nano-corr {input} --out-dir assembly/ --genome-size 0.2m --meta -t 8"

rule medaka:
    input:
        "data/samples/{sample}.fastq"
    output:
        "medaka_output/{sample}_medaka.fa"
    conda:
        "envs/medaka.yaml"
    shell:
        "medaka_consensus -f -i {input} -d {input.prev_fa} -o assemblies/{wildcards.sample}_{wildcards.assembly}+medaka -t {threads}"


rule mapping_reads:
    input:
        "trimmed/{sample}.trimmed.fastq"
    output:
        temp("mapped/{sample}.sam")
    conda:
        "envs/mapping.yaml"
    shell:
        "minimap2 -ax map-ont data/monkeypox.fa {input} > {output}"

rule sam_to_bam:
    input:
        "mapped/{sample}.sam"
    output:
        temp("mapped/{sample}.bam")
    conda:
        "envs/mapping.yaml"
    shell:
        "samtools view -S -b {input} > {output}"

rule sort:
    input:
        "mapped/{sample}.bam"
    output:
        "mapped/{sample}.sorted.bam"
    conda:
        "envs/mapping.yaml"
    shell: 
        "samtools sort {input} -o {output}"

rule index_bam:
    input:
        "mapped/{sample}.sorted.bam"
    output:
        "mapped/{sample}.sorted.bam.bai"
    conda:
        "envs/mapping.yaml"
    shell:
        "samtools index {input} {output}"

rule consensus_generation:
    input:
        "mapped/{sample}.sorted.bam"
    output:
        "consensus/{sample}.fa"
    conda:
        "envs/mapping.yaml"
    shell:
        "samtools mpileup -aa -A -Q 10 {input} | ivar consensus -q 8 -m 8 -p {output}"


rule bcftools_call:
    input:
        fa="data/monkeypox.fa",
        reads="mapped/{sample}.sorted.bam"
    output:
        "calls/{sample}.vcf"
    conda:
        "envs/calling.yaml"
    shell:
        "freebayes -f {input.fa} -C 5 -F 0.01 -b {input.reads} >{output}"

rule plot_qual:
    input:
        "calls/{sample}.vcf"
    output:
        plot_pdf="plots/quals.{sample}.pdf"
    conda:
        "envs/plot-quals.yaml"
    script:
        "scripts/plot-quals.py"


rule homopolish:
    conda: "env/conda-homopolish.yaml"
    input:
        prev_fa = "consensus/{sample}.fa",
        ref = "references/{ref}.fa"
    output: 
        fa = "assemblies/{sample}_{assembly}+homopolish/output_{ref}.fa",
    shell:
        "homopolish polish -a {input.prev_fa} -m R9.4.pkl -o polished/{}.fa -l {input.ref}"

rule busco:
    input:
        "consensus/{sample}.fa"
    output:
        "busco_output/{sample}_busco.txt"
    conda:
        "envs/busco.yaml"
    shell:
        "busco -f -c 20 -m genome -i input -o {output} --auto-lineage-prok"
