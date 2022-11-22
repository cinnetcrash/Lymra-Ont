rule map_reads:
    input:
        "data/genome.fa"
    output:
        "mapped/{output}.bam"
    conda:
        "envs/mapping.yaml"
    shell:
        "bwa mem {input} | samtools view -Sb - > {output}"

rule sort:
    input:
        "mapped/{sample}.bam"
    output:
        ""