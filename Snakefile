configfile: "config.yaml"


rule all:
    "busco_output/{sample}_busco.txt"

rule porechop_trim:
    input:
        "data/samples/{sample}.fastq.gz"
    output:
        "trimmed/{sample}.trimmed.fastq"
    conda: 
        "envs/conda-porechop.yaml"
    shell: 
        "porechop -i {input} -o {output} -t 4 --barcode_threshold 60 --barcode_diff 1"

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
        "kraken2 --db {input.database} {input.trimmed_reads}--report {output.report_kraken2} --classified-out {output.classified_out}"

rule flye:
    input:
        viral_reads="trimmed/{sample}_viral_reads.fastq",
        output_dir="assembly/flye/"
    output:
        "assembly/{sample}.fasta"
    conda:
        "envs/conda-flye.yaml"
    shell:
        "flye --nano-corr {input.viral_reads} --out-dir {output.output_dir} --genome-size 0.2m --meta -t 8"

rule medaka:
    input:
        fq="trimmed/{sample}_viral_reads.fastq",
        reference="data/monkeypox.fa"
    output:
        "medaka_output/{sample}_medaka.fasta"
    conda:
        "envs/conda-medaka.yaml"
    shell:
        "medaka_consensus -i {input.fq} -d {input.reference} -o {output} -t 8 -m r941_min_high_g303"

rule homopolish:
    conda: "envs/conda-homopolish.yaml"
    input:
        prev_fa = "medaka_output/{sample}_medaka.fasta",
        ref = "data/monkeypox.fa"
    output: 
        fa = "polish/{sample}_homopolish/{sample}_polished.fa"
    shell:
        "homopolish polish -a {input.prev_fa} -m R9.4.pkl -o {output.fa} -l {input.ref}"

rule busco:
    input:
        "polish/{sample}_homopolish/{sample}_polished.fa"
    output:
        "busco_output/{sample}_busco.txt"
    conda:
        "envs/conda-porechop.yaml"
    shell:
        "busco -f -c 20 -m genome -i {input} -o {output} --auto-lineage-prok"
