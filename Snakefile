import os
from snakemake.io import expand

# This snakefile contains the processing steps of analysing 
# the 2x75 bp paired end of Lagomarsino et al. 2021. 
# 
# define some variables. 
pair_ids = [".r1", ".r2"]
ext = ".fq.gz"
FQC_EXT = ["zip", "html"]

# some python code to extract the "sample id" from the file names. 
# Lachie wrote this for me since I'm not a snake queen yet 
samples = os.listdir("01_rawdata/fastqraw")
samples = [sample.replace(ext, "") for sample in samples]
for id in pair_ids:
    samples = [sample.replace(id, "") for sample in samples]
samples = list(set(samples))

# here, we need to define the files which are the end products of the workflow. 
# Any files which are the input for a next step do not need to be included. 
# I have all of them in here, since I built this pipeline from scratch essentially. 
# This expand function will generate file names which alter based on these "wildcardds". 
# A nice explanantion of how expand() works can be found here:
# https://endrebak.gitbooks.io/the-snakemake-book/content/chapters/expand/expand.html 

rule all:
    input:
        expand("01_rawdata/fastqc/{SAMPLE}{PAIR}_fastqc.{EXT}", SAMPLE = samples, PAIR = pair_ids, EXT = FQC_EXT),
        expand("02_trimdata/fastq/{SAMPLE}{PAIR_ID}{EXT}", SAMPLE = samples, PAIR_ID = pair_ids, EXT = ext),
        expand("02_trimdata/fastqc/{SAMPLE}{PAIR}_fastqc.{EXT}", SAMPLE = samples, PAIR = pair_ids, EXT = FQC_EXT),
        expand("03_alignstar/bam/{SAMPLE}.Aligned.sortedByCoord.out.bam", SAMPLE = samples),
        expand("03_alignstar/bam/{SAMPLE}.Aligned.sortedByCoord.out.bam.bai", SAMPLE = samples),
        expand("03_alignstar/FastQC/{SAMPLE}.Aligned.sortedByCoord.out_fastqc.{EXT}", SAMPLE = samples, EXT = FQC_EXT), 
        "04_featureCounts/counts.out"

# First, I will run fastqc on the raw data. 
rule fastqc:
    input:
        R1 = "01_rawdata/fastqraw/{SAMPLE}.fq.gz"
    params:
        outdir = "01_rawdata/fastqc/"
    output:    
        html = "01_rawdata/fastqc/{SAMPLE}_fastqc.html",
        zip = "01_rawdata/fastqc/{SAMPLE}_fastqc.zip"
    conda:
        "smk/envs/default.yaml"
    resources: # parameters which will submit to phoenix
        cpu = 1,
        ntasks = 1,
        time = "00-01:00:00",
        mem_mb = 4000
    shell:
        """
        fastqc \
        -t {resources.cpu} \
        -o {params.outdir} \
        {input}
        """

# The next step is to run adaptor removal using fastp. 
# I will only retain reads which are more than 20 nt in length after 
# trimming, and have a quality score of at least 15 phred. 
# I will also discard polyG reads, which I noticed are in this dataset. 

rule trim:
    input:
        R1 = "01_rawdata/fastqraw/{SAMPLE}.r1.fq.gz",
        R2 = "01_rawdata/fastqraw/{SAMPLE}.r2.fq.gz"
    output:
        R1 = "02_trimdata/fastq/{SAMPLE}.r1.fq.gz",
        R2 = "02_trimdata/fastq/{SAMPLE}.r2.fq.gz",
        json = "02_trimdata/log/{SAMPLE}.json",
        html = "02_trimdata/log/{SAMPLE}.html"
    params:
        bname = "02_trimdata/fastq/{SAMPLE}"
    conda:
        "smk/envs/default.yaml"
    resources:
        cpu = 1,
        ntasks = 1,
        time = "00-01:00:00",
        mem_mb = 4000
    shell:
        """        
        fastp \
            -l 20 \
            --json {output.json} \
            --html {output.html} \
            --detect_adapter_for_pe \
            --trim_poly_g \
            --out1 {output.R1} \
            --out2 {output.R2} \
            -i {input.R1} \
            -I {input.R2}
        """
    # input:
    #     R1 = "01_rawdata/fastqraw/{SAMPLE}.r1.fq.gz",
    #     R2 = "01_rawdata/fastqraw/{SAMPLE}.r2.fq.gz"
    # output:
    #     R1 = "02_trimdata/fastq/{SAMPLE}.r1.fq.gz",
    #     R2 = "02_trimdata/fastq/{SAMPLE}.r2.fq.gz"
    # params:
    #     bname = "02_trimdata/fastq/{SAMPLE}"
    # conda:
    #     "smk/envs/default.yaml"
    # resources:
    #     cpu = 1,
    #     ntasks = 1,
    #     time = "00-01:00:00",
    #     mem_mb = 4000
    # shell:
    #     """        
    #     AdapterRemoval \
    #         --gzip \
    #         --trimns \
    #         --trimqualities \
    #         --minquality 30 \
    #         --minlength 20 \
    #         --basename {params.bname} \
    #         --threads {resources.cpu} \
    #         --output1 {output.R1} \
    #         --output2 {output.R2} \
    #         --file1 {input.R1} \
    #         --file2 {input.R2}
    #     """
# repeat fastqc after trimming.         
rule trimqc:
    input:
        R1 = "02_trimdata/fastq/{SAMPLE}.fq.gz"
    params:
        outdir = "02_trimdata/fastqc/"  
    output:    
        html = "02_trimdata/fastqc/{SAMPLE}_fastqc.html",
        zip = "02_trimdata/fastqc/{SAMPLE}_fastqc.zip"
    conda:
        "smk/envs/default.yaml"
    resources:
        cpu = 2,
        ntasks = 2,
        time = "00-01:00:00",
        mem_mb = 6000
    shell:
        """
        fastqc \
        -f fastq \
        -t {resources.cpu} \
        -o {params.outdir} \
        --noextract \
        {input}
        """
        
rule align:
# here, we are aligning to the human genome (ensembll release 98). 
# The genome index is already generated in the /hpcfs folder by 
    input:
        R1 = "02_trimdata/fastq/{SAMPLE}.r1.fq.gz", 
        R2 = "02_trimdata/fastq/{SAMPLE}.r2.fq.gz"       
    params:
        genomedir = "/hpcfs/archive/biorefs/reference_genomes/ensembl-release-98/homo_sapiens/star/",
        bname = "03_alignstar/bam/{SAMPLE}."
    output:    
        bam = "03_alignstar/bam/{SAMPLE}.Aligned.sortedByCoord.out.bam"
    conda:
        "smk/envs/default.yaml"
    resources:
        cpu = 16,
        ntasks = 1,
        time = "00-05:00:00",
        mem_mb = 50000
    shell:
        """
        STAR \
        --genomeDir {params.genomedir}\
        --runThreadN {resources.cpu} \
        --readFilesIn {input.R1} {input.R2} \
        --readFilesCommand "gunzip -c" \
        --outSAMtype BAM SortedByCoordinate \
        --outFileNamePrefix {params.bname}

        mkdir -p 03_alignstar/log
		mv {params.bname}*out 03_alignstar/log
		mv {params.bname}*tab 03_alignstar/log
        """

rule fastqcalign:
    input:
        bam = "03_alignstar/bam/{SAMPLE}.Aligned.sortedByCoord.out.bam"
    output:
        "03_alignstar/FastQC/{SAMPLE}.Aligned.sortedByCoord.out_fastqc.zip",
        "03_alignstar/FastQC/{SAMPLE}.Aligned.sortedByCoord.out_fastqc.html"
    params:
        outDir = "03_alignstar/FastQC/"
    conda:
        "smk/envs/default.yaml"
    resources:
        cpu = 2,
        ntasks = 2,
        time = "00-01:00:00",
        mem_mb = 6000
    shell:
        """
        fastqc -t {resources.cpu} \
        -o {params.outDir} \
        --noextract \
        {input.bam}
        """

rule indexBam:
    input:
        "03_alignstar/bam/{SAMPLE}.Aligned.sortedByCoord.out.bam"
    output:
        "03_alignstar/bam/{SAMPLE}.Aligned.sortedByCoord.out.bam.bai"
    conda:
        "smk/envs/default.yaml"
    resources:
        cpu = 2,
        ntasks = 1, 
        time = "00-01:00:00",
        mem_mb = 6000
    shell:
        """
        samtools index {input} {output}
        """


rule featureCounts:
    input:
        bam = expand("03_alignstar/bam/{SAMPLE}.Aligned.sortedByCoord.out.bam", SAMPLE = samples),
        gtf = "Homo_sapiens.GRCh38.98.chr.gtf" 
    output:
        counts = "04_featureCounts/counts.out"
    conda:
        "smk/envs/default.yaml"
    resources:
        cpu = 2,
        ntasks = 2,
        time = "00-01:00:00",
        mem_mb = 6000
    params:
        minOverlap = 1,
        fracOverlap = 1,
        q = 10
# Settings for featureCounts. 
# To set this to count strictly exonic reads, I change fracOverlap to be the value 1. 
# The value minOverlap may also need adjusting based on your own read lengths. 
# the -Q option is set to 10, meaining a mapping quality of at least 10. 
# the -p flad indicates the input bam files were generated from paired end data
    shell:
       """
       featureCounts \
       -Q {params.q} \
       -p \
       --minOverlap {params.minOverlap} \
       --fracOverlap {params.fracOverlap} \
       -T {resources.cpu} \
       -a {input.gtf} \
       -o {output} \
       {input.bam}
       """
