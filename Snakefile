SAMPLES = ["Leaf51", "Axenic", "Leaf49", "Leaf137", "Leaf405", "Leaf177", "Leaf257"]
EXPANSION = ["1", "2", "3", "4", "5"]
INDICES = ["TAIR10", "ARA11"]
THREADS = 8

from Scripts.get_samples import get_sample_entries

entries = get_sample_entries("Scripts/samples.csv")
ALL_SAMPLES = {}
for name, entry in entries.items():
    if name in SAMPLES:
        ALL_SAMPLES[name] = [sample.submitted_ftp.split("/")[-1] for sample in entry]


rule all:
    input:
        expand("output/hisat2/{sample}_{expansion}.bam.bai", sample=SAMPLES, expansion=EXPANSION, indices=INDICES),
        expand("output/prepDE/{indices}_gene_count_matrix.csv", indices=INDICES),


rule hisat2:
    input:
        reads= lambda wildcards: expand("Samples/transcriptome/{file}", file=ALL_SAMPLES[wildcards.sample][int(wildcards.expansion)-1]),
        index="Samples/genomes/Arabidopsis_thaliana/TAIR.1.ht2" # this file is just there to make sure the index is present
    output:
        temp("output/hisat2/{sample}_{expansion}.sam") 
    params:
        index="Samples/genomes/Arabidopsis_thaliana/TAIR" 
    threads: THREADS # to avoid too much load on Altschul
    shell:
        "hisat2 -p {threads} -x {params.index} -U {input.reads} -S {output} --dta"


rule sort_sam:
    input:
        "output/hisat2/{sample}_{expansion}.sam"
    output:
        "output/hisat2/{sample}_{expansion}.bam"
    threads: THREADS
    shell:
        "samtools sort -o {output} {input} -@ {threads}"

rule index_bam:
    input:
        "output/hisat2/{sample}_{expansion}.bam"
    output:
        "output/hisat2/{sample}_{expansion}.bam.bai"
    shell:
        "samtools index {input}"

rule stringtie_TAIR10:
    input:
        bamfile="output/hisat2/{sample}_{expansion}.bam", # stringtie is only used after hisat2
        annotation="Samples/genomes/Arabidopsis_thaliana/genes.gtf"
    output:
        "output/stringtie/TAIR10/{sample}_{expansion}/{sample}_{expansion}.gtf"
    params:
        label="{sample}_{expansion}"
    shell:
        "stringtie -e -G {input.annotation} -o {output} -l {params.label} {input.bamfile}"

rule stringtie_ARA11:
    input:
        bamfile="output/hisat2/{sample}_{expansion}.bam", # stringtie is only used after hisat2
        annotation="Samples/genomes/Arabidopsis_thaliana/Araport11_GTF_cleaned_chromosome_names.gtf"
    output:
        "output/stringtie/ARA11/{sample}_{expansion}/{sample}_{expansion}.gtf"
    params:
        label="{sample}_{expansion}"
    threads: THREADS
    shell:
        "stringtie -e -G {input.annotation} -o {output} -l {params.label} {input.bamfile} -p {threads}"

rule prepDE:
    input:
        expand("output/stringtie/{indices}/{sample}_{expansion}/{sample}_{expansion}.gtf", sample=SAMPLES, expansion=EXPANSION, indices=INDICES)
    output:
        "output/prepDE/{indices}_gene_count_matrix.csv",
        "output/prepDE/{indices}_transcript_count_matrix.csv"
    shell:
        "python3 Scripts/prepDE.py -i output/stringtie/{wildcards.indices} -g output/prepDE/{wildcards.indices}_gene_count_matrix.csv -t output/prepDE/{wildcards.indices}_transcript_count_matrix.csv"



