#!/usr/bin/env python
from os.path import join, abspath, expanduser
import re

# Read in config file
OUTDIR = config['out_dir']
READS_DIR = config['reads_dir']
SAMP = config['samp_prefixes']


# Convert out directory to absolute path
if OUTDIR[0] == '~':
    OUTDIR = expanduser(OUTDIR)
OUTDIR = abspath(OUTDIR)

# Get fastq file names
FILES = [f for f in os.listdir(READS_DIR) if f.endswith(".fastq")]
SAMPLE_PREFIX = list(re.split("_[12].fastq.gz", i)[0] for i in FILES)


###############################################################

rule all:
    input:
        join(OUTDIR, "RNA_preprocess/00_pre_fastqc/multiqc_report.html"),
        join(OUTDIR, "RNA_preprocess/02_post_fastqc/multiqc_report.html"),

###############################################################
################### Quality check and trim ####################
###############################################################

# quality check raw reads
rule pre_fastqc:
    input:
        reads = join(READS_DIR, "{sample}_R{R}_001.fastq.gz")
    output:
        join(OUTDIR, "RNA_preprocess/00_pre_fastqc/{sample}_R{R}_001_fastqc.html")
    params:
        outdir = join(OUTDIR, "RNA_preprocess/00_pre_fastqc/")
    shell: """
        module load fastqc/0.12.1

        mkdir -p {params.outdir}
        fastqc {input} --outdir {params.outdir}
    """

# compile fastqc files into single report
rule pre_multiqc:
    input:
        expand(join(OUTDIR, "RNA_preprocess/00_pre_fastqc/{sample}_R{R}_001_fastqc.html"), sample=SAMP, R=["1","2"])
    output:
        join(OUTDIR, "RNA_preprocess/00_pre_fastqc/multiqc_report.html")
    params:
        dir = join(OUTDIR, "RNA_preprocess/00_pre_fastqc/")
    shell:"""
        module load multiqc/1.25
        multiqc --force {params.dir} -o {params.dir}
    """

# Trim reads
rule fastp:
    input:
        r1 = join(READS_DIR, "{sample}_R1_001.fastq.gz"),
        r2 = join(READS_DIR, "{sample}_R2_001.fastq.gz")
    output:
        r1 = join(OUTDIR, "RNA_preprocess/01_fastp_trim/{sample}_R1.trimmed.fastq.gz"),
        r2 = join(OUTDIR, "RNA_preprocess/01_fastp_trim/{sample}_R2.trimmed.fastq.gz"),
        html = join(OUTDIR, "RNA_preprocess/01_fastp_trim/{sample}_fastp.html"),
        json = join(OUTDIR, "RNA_preprocess/01_fastp_trim/{sample}_fastp.json"),
    resources:
        mem = 24,
        time = 4
    threads: 8
    shell: """
        module load fastp/0.23.4
        # Trim first 13 bp of all reads and default adapter trimming
        fastp --trim_front1 13 --trim_front2 13 -i {input.r1} -I {input.r2} -o {output.r1} -O {output.r2} -h {output.html} -j {output.json} -w {threads}
    """

# quality check raw reads
rule post_fastqc:
    input:
        reads = join(OUTDIR, "RNA_preprocess/01_fastp_trim/{sample}_{R}.trimmed.fastq.gz")
    output:
        join(OUTDIR, "RNA_preprocess/02_post_fastqc/{sample}_{R}.trimmed_fastqc.html")
    params:
        outdir = join(OUTDIR, "RNA_preprocess/02_post_fastqc/")
    shell: """
        module load fastqc/0.12.1

        mkdir -p {params.outdir}
        fastqc {input} --outdir {params.outdir}
    """

# compile fastqc files into single report
rule post_multiqc:
    input:
        expand(join(OUTDIR, "RNA_preprocess/02_post_fastqc/{sample}_R{R}.trimmed_fastqc.html"), sample=SAMP, R=["1","2"])
    output:
        join(OUTDIR, "RNA_preprocess/02_post_fastqc/multiqc_report.html")
    params:
        dir = join(OUTDIR, "RNA_preprocess/02_post_fastqc/")
    shell:"""
        module load multiqc/1.25
        multiqc --force {params.dir} -o {params.dir}
    """