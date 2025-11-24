import subprocess
import os
import pandas as pd

configfile: "config/config.yaml"

samplefile = 'copied_reads_dir/samples_to_trim.csv'
samples_df = pd.read_csv(samplefile,sep = '\t')

SAMPLE = samples_df['Sample'].tolist()

rule all:
    input:
        r1 = expand("copied_reads_dir/trimmed_reads/{sample}_R1_trim_paired.fastq.gz",sample=SAMPLE),
        r2 = expand("copied_reads_dir/trimmed_reads/{sample}_R2_trim_paired.fastq.gz",sample=SAMPLE),

rule trimmomatic_pe:
    input:
        r1 = "copied_reads_dir/raw_reads/{sample}_R1.fastq.gz",
        r2 = "copied_reads_dir/raw_reads/{sample}_R2.fastq.gz",
    output:
        r1 = "copied_reads_dir/trimmed_reads/{sample}_R1_trim_paired.fastq.gz",
        r2 = "copied_reads_dir/trimmed_reads/{sample}_R2_trim_paired.fastq.gz", 
        # reads where trimming entirely removed the mate
        r1_unpaired = "copied_reads_dir/temp/{sample}/{sample}_R1_trim_unpaired.fastq.gz",
        r2_unpaired = "copied_reads_dir/temp/{sample}/{sample}_R2_trim_unpaired.fastq.gz",
    params:
        adapter_filepath=config["adapter_file"],
        seed=config["seed_mismatches"],
        palindrome_clip=config["palindrome_clipthreshold"],
        simple_clip=config["simple_clipthreshold"],
        minadapterlength=config["minadapterlength"],
        keep_both_reads=config["keep_both_reads"],
        window_size=config["window_size"],
        window_size_quality=config["window_size_quality"],
        minlength=config["minlength"],
        headcrop_length=config["headcrop_length"],
        lead_trail_qual=config["lead_trail_qual"],
    log:
        "copied_reads_dir/temp/{sample}/{sample}.log"
    threads: 4
    singularity:
        "docker://staphb/trimmomatic:0.39"
    #retries: 1
    resources:
        mem_mb = 15000,
        #mem_mb = lambda wildcards, attempt: 15000 + ((attempt-1)*15000),
        runtime = 120,
    shell:
        """
        trimmomatic PE {input.r1} {input.r2} {output.r1} {output.r1_unpaired} {output.r2} {output.r2_unpaired} -threads {threads} \
        ILLUMINACLIP:{params.adapter_filepath}:{params.seed}:{params.palindrome_clip}:{params.simple_clip}:{params.minadapterlength}:{params.keep_both_reads} \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:{params.window_size}:{params.window_size_quality} MINLEN:{params.minlength} HEADCROP:{params.headcrop_length} &>{log}
        """
