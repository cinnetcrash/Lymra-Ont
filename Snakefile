rule map_reads:
    input:
        "data/monkeypox.fa",
        "../{sample}.fastq"
    output:
        "mapped/{sample}.bam"
    conda:
        "envs/mapping.yaml"
    shell:
        "bwa mem {input} | samtools view -Sb - > {output}"

rule sort:
    input:
        "mapped/{sample}.bam"
    output:
        "mapped/{sample}.sorted.bam"
    conda:
        "envs/mapping.yaml"
    shell: 
        "samtools sort -o {output} {input}"