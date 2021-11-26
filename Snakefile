SAMPLES = ["X"]
TOOLS = ["bowtie2", "hisat2"]

rule all:
    input:
        expand("{sample}_{tool}.bam.bai", sample=SAMPLES, tool=TOOLS),
        expand("{sample}_hisat2.gtf", sample=SAMPLES),
        expand("{sample}_kallisto/abundance.tsv", sample=SAMPLES),
        expand("{sample}_trinity.blast.gff", sample=SAMPLES)

rule bowtie2:
    input:
        reads="{sample}.fastq",
        index="Bowtie2Index/genome.1.bt2" # this file is just there to make sure the index is present
    output:
        temp("{sample}_bowtie2.sam") # temp() removes the sam file when the pipeline is done
    params:
        index="Bowtie2Index/genome" # the index is actually a group of files with this prefix
    shell:
        "bowtie2 -x {params.index} -r {input.reads} -S {output}"

rule hisat2:
    input:
        reads="{sample}.fastq",
        index="TAIR.1.ht2" # this file is just there to make sure the index is present
    output:
        temp("{sample}_hisat2.sam") # temp() removes the sam file when the pipeline is done
    params:
        index="TAIR" # the index is actually a group of files with this prefix
    threads: 1 # to avoid too much load on Altschul
    shell:
        "hisat2 -p {threads} -x {params.index} -q {input.reads} -S {output}"

rule sort_sam: # works for both bowtie2 and hisat2
    input:
        "{sample}_{tool}.sam"
    output:
        "{sample}_{tool}.bam"
    shell:
        "samtools sort -o {output} {input}"

rule index_bam: # works for both bowtie2 and hisat2
    input:
        "{sample}_{tool}.bam"
    output:
        "{sample}_{tool}.bam.bai"
    shell:
        "samtools index {input}"

rule stringtie:
    input:
        bamfile="{sample}_hisat2.bam", # stringtie is only used after hisat2
        annotation="genes.gtf"
    output:
        "{sample}_hisat2.gtf"
    params:
        label="{sample}"
    shell:
        "stringtie -G {input.annotation} -o {output} -l {params.label} {input.bamfile}"

rule trinity:
    input:
        "{sample}.fastq"
    output:
        "{sample}_trinity/Trinity.fasta"
    params:
        seqtype="fq",
        maxmem="10G",
        outdir="{sample}_trinity"
    shell:
        "Trinity --seqType {params.seqtype} --max_memory {params.maxmem} --output {params.outdir} --single {input}"

rule kallisto:
    input:
        reads="{sample}.fastq",
        index="TAIR.idx"
    output:
        "{sample}_kallisto/abundance.tsv"
    params:
        length=200,
        sd=20,
        outdir="{sample}_kallisto"
    threads: 1
    shell:
        "kallisto quant -i {input.index} -o {params.outdir} --single -l {params.length} -s {params.sd} -t {threads} {input.reads}"

rule make_blastdb:
    input:
        "genome.fa"
    output:
        "TAIR10.nhr" # just one of the files of the blast database
    params:
        db="TAIR10",
        type="nucl"
    shell:
        "makeblastdb -in {input} -dbtype {params.type} -out {params.db}"

rule blastn:
    input:
        query="{sample}_trinity/Trinity.fasta/Trinity.fasta",
        db="TAIR10.nhr"
    output:
        file="{sample}_trinity.blast"
    params:
        db="TAIR10",
        evalue="1e-10",
        format=7
    shell:
        "blastn -query {input.query} -db {params.db} -evalue {params.evalue} -outfmt {params.format} -out {output.file}"

rule blast2gff:
    input:
        "{sample}_trinity.blast"
    output:
        "{sample}_trinity.blast.gff"
    shell:
        "Rscript blast2gff.R {input}"

