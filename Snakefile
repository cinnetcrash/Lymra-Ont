configfile: "config.yaml"

rule all:
    "busco_output/"

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
        "kraken2 --db {input.database} {input.trimmed_reads} --report {output.report_kraken2} --classified-out {output.classified_out}"

rule flye:
    input:
        viral_reads="trimmed/{sample}_viral_reads.fastq",
    output:
        "assembly/{sample}.fasta"
    params:
        output_dir="assembly/flye/"
    conda: srcdir("envs/conda-flye.yaml")
    shell:
        "flye --nano-corr {input.viral_reads} --out-dir {params.output_dir} --genome-size 0.2m --meta -t 8"

rule minimap2_sort:
    input:
        fastq="data/samples/{sample}.fastq.gz",
        reference=config["reference"]
    params:
        preset="-x map-ont", 
        kmer="-k 21", 
        secondary_aligment="--secondary=yes -N 5", 
        sec_score = "-p 0.8" 
    output:
        bam="data/mapped/{sample}.sorted.bam",
        bai="data/mapped/{sample}.sorted.bam.bai"
    conda: srcdir("envs/conda-flye.yaml")
    shell:
        """
        minimap2 -a --frag=yes -t 6 {params.preset} {params.kmer} {params.secondary_aligment} \
        {params.sec_score} {input.reference} {input.fastq} | samtools view -bS - | samtools sort -o {output.bam} - ; \n
        samtools index {output.bam}
        """

rule medaka1:
    input:
        bam="data/mapped/{sample}.sorted.bam",
        reference=config["reference"]
    params:
        model=config["medaka_model"]
    conda: srcdir("envs/conda-medaka.yaml")
    output:
        hdf="data/medaka/{sample}.sorted.hdf"
    shell:
        "medaka consensus {input.bam} {output.hdf}"
        
rule medaka2:
    input:
        hdf_file="data/medaka/{sample}.sorted.hdf",
        reference=config["reference"]
    params:
        model=config["medaka_model"]
    conda: srcdir("envs/conda-medaka.yaml")
    output:
        vcf="data/medaka/{sample}.sorted.vcf"
    shell:
       "medaka variant {input.reference} {input.hdf_file} {output.vcf}"

rule tabix_medaka:
    input:
        vcf="data/medaka/{sample}.sorted.vcf"
    conda: srcdir("envs/conda-porechop.yaml")
    output:
        gz="data/medaka/{sample}.sorted.vcf.gz"
    shell:
        "gzip -f {input.vcf} | tabix -p vcf {output.gz}"

rule longshot:
    input:
        bam="data/mapped/{sample}.sorted.bam",
        vcf="data/medaka/{sample}.sorted.vcf",
        reference=config["reference"]
    params:
        pvalue="-P 0", 
        flags="-F -A -n "
    conda: srcdir("envs/conda-medaka.yaml")
    output:
        vcf="data/longshot/{sample}.vcf"
    shell:
        "medaka tools annotate --pad 25 {input.vcf} {input.reference} {input.bam} {output.vcf}"

rule qualimap:
    input:
        bam="data/mapped/{sample}.sorted.bam"
    params:
        java="--java-mem-size=16G"
    conda: srcdir("envs/conda-porechop.yaml")
    output:
        html=directory("data/qualimap/{sample}")
    shell:
        "qualimap bamqc -bam {input.bam} {params.java} --outdir {output} --outformat html"

rule margin_cons_medaka:
    input:
        vcf="data/longshot/{sample}.vcf",
        bam="data/mapped/{sample}.sorted.bam",
        reference=config["reference"]
    params:
        depth=config["coverage"],
        qual=config["quality"]
    conda: srcdir("envs/ncov.yml")
    output:
        fasta="data/genome/{sample}.consensus.fasta",
        report="data/genome/{sample}.report.txt"
    script:
        "scripts/margin_cons_medaka.py"

rule busco:
    input:
        baam="data/mapped/{sample}.sorted.bam"
    output:
        "data"
    params:
        exit="busco_output/"
    conda:
        "envs/conda-porechop.yaml"
    shell:
        "busco -i {input.baam} -o {params.exit} --auto-lineage-prok -m genome"

