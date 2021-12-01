SAMPLES = ["Leaf51", "Axenic", "Leaf49", "Leaf137", "Leaf405", "Leaf177", "Leaf257"]
EXPANSION = ["1", "2", "3"]

rule all:
    input:
        expand("output/hisat2/{sample}_{expansion}.bam.bai", sample=SAMPLES, expansion=EXPANSION),
        ".lock"


rule hisat2:
    input:
        reads="Samples/{sample}_{expansion}.fastq.gz",
        index="Samples/genomes/Arabidopsis_thaliana/TAIR.1.ht2" # this file is just there to make sure the index is present
    output:
        temp("output/hisat2/{sample}_{expansion}.sam") \
        # temp() removes the sam file when the pipeline is done
    params:
        index="Samples/genomes/Arabidopsis_thaliana/TAIR" \
        # the index is actually a group of files with this prefix
    threads: 1 # to avoid too much load on Altschul
    shell:
        "hisat2 -p {threads} -x {params.index} -U {input.reads} -S {output} --dta"

rule sort_sam:
    input:
        "output/hisat2/{sample}_{expansion}.sam"
    output:
        "output/hisat2/{sample}_{expansion}.bam"
    shell:
        "samtools sort -o {output} {input}"

rule index_bam:
    input:
        "output/hisat2/{sample}_{expansion}.bam"
    output:
        "output/hisat2/{sample}_{expansion}.bam.bai"
    shell:
        "samtools index {input}"

rule stringtie:
    input:
        bamfile="output/hisat2/{sample}_{expansion}.bam", # stringtie is only used after hisat2
        #annotation="Samples/genomes/Arabidopsis_thaliana/Araport11_GTF_cleaned_chromosome_names.gtf"
        annotation="Samples/genomes/Arabidopsis_thaliana/genes.gtf"
    output:
        "output/stringtie/{sample}_{expansion}/{sample}_{expansion}.gtf"
    params:
        label="{sample}_{expansion}"
    shell:
        "stringtie -e -G {input.annotation} -o {output} -l {params.label} {input.bamfile}"

rule prepDE:
    input:
        expand("output/stringtie/{sample}_{expansion}/{sample}_{expansion}.gtf", sample=SAMPLES, expansion=EXPANSION)
    output:
        "output/prepDE/gene_count_matrix.csv",
        "output/prepDE/transcript_count_matrix.csv"
    shell:
        "python3 Scripts/prepDE.py -i output/stringtie -g output/prepDE/gene_count_matrix.csv -t output/prepDE/transcript_count_matrix.csv"

rule DESeq:
    input:
        "output/prepDE/gene_count_matrix.csv"
    output:
        ".lock"
    shell:
        "Rscript"