import subprocess
import pandas as pd
import os

configfile: "config/config.yaml"

samples_df = pd.read_csv(config["samples"])
SAMPLE = list(samples_df['sample_id'])

PREFIX = config["prefix"]

# prepare reference directory
# note that difference reference genomes are not supported yet - the reference must be the same for all samples!
# this directory also must contain a .gff file for the reference genome
def prep_reference(ref_fasta,ref_gff):
    ref_name = ref_fasta.split('/')[-1].split('.fa')[0]
    ref_dir = f'results/reference/{ref_name}/'
    new_ref_fasta = f'results/reference/{ref_name}/{ref_name}.fasta'
    new_ref_gff = f'results/reference/{ref_name}/{ref_name}.gff'
    print(f'The reference used for ALL samples in this run is {ref_name}, located at {ref_fasta}')
    if not os.path.exists(new_ref_fasta) or not os.path.exists(new_ref_gff):
        subprocess.call(['mkdir','-p',ref_dir])
        subprocess.call(['cp',ref_fasta,new_ref_fasta])
        subprocess.call(['cp',ref_gff,new_ref_gff])
        print(f'A new reference directory has been created at {ref_dir}')
    else:
        print('Using existing reference files at this location!')
    print('Remember to check the reference genome mitochondrial chromosome location in config/config.yaml')
    return([new_ref_fasta,new_ref_gff,ref_dir])

# this will copy the provided reference to a new directory (assuming it's not already there)
REFERENCE,REF_GFF,REF_DIR = prep_reference(config["reference_fasta"],config["reference_gff"])

rule all:
    input:
        annotated_variants = expand("results/{prefix}/{sample}/annotate_SVs/annotated_variants.tab",sample=SAMPLE,prefix=PREFIX),
        #call_svs = expand("results/{prefix}/{sample}/call_SVs/perSVade_finished_vaf_file.txt",sample=SAMPLE,prefix=PREFIX),

rule ref_infer_repeats:
    input:
        ref = REFERENCE
    output:
        repeats = REF_DIR + 'repeat_inference/combined_repeats.tab',
    params:
        outdir = REF_DIR + 'repeat_inference/',
        fraction_available_mem = lambda wildcards, resources: round((resources.mem_mb/1000)/175 - 0.005,2),
    resources:
        mem_mb=10000,
        runtime=240,
    threads: 4
    singularity:
        "docker://mikischikora/persvade:v1.02.6"
    shell:
        """
        set +u
        source /opt/conda/etc/profile.d/conda.sh
        conda activate perSVade_env
        python /perSVade/scripts/perSVade infer_repeats --ref {input.ref} -o {params.outdir} \
        --replace --verbose --fraction_available_mem {params.fraction_available_mem} --threads {threads}
        set -u
        """


# these steps have been changed to remove the trimmed reads once alignment is done

rule step1_qc:
    input:
        r1 = config["short_reads"] + "/" + "{sample}_R1_trim_paired.fastq.gz",
        r2 = config["short_reads"] + "/" + "{sample}_R2_trim_paired.fastq.gz",
    output:
        #r1trim = "results/{sample}/misc/trimmed_reads/trimmed_reads1.fastq.gz",
        #r2trim = "results/{sample}/misc/trimmed_reads/trimmed_reads2.fastq.gz",
        trim_finish = "results/{prefix}/{sample}/misc/trimmed_reads/perSVade_finished_file.txt",
    params:
        outdir = "results/{prefix}/{sample}/misc/trimmed_reads/",
        fraction_available_mem = lambda wildcards, resources: round((resources.mem_mb/1000)/175 - 0.005,2),
    resources:
        mem_mb=5000,
        runtime=360,
    threads: 1
    singularity:
        "docker://mikischikora/persvade:v1.02.6"
    shell:
        """
        set +u
        source /opt/conda/etc/profile.d/conda.sh
        conda activate perSVade_env
        python /perSVade/scripts/perSVade trim_reads_and_QC -f1 {input.r1} -f2 {input.r2} -o {params.outdir} \
        --replace --verbose --fraction_available_mem {params.fraction_available_mem} --threads {threads}
        set -u
        """

rule step2_align:
    input:
        #r1trim = "results/{sample}/misc/trimmed_reads/trimmed_reads1.fastq.gz",
        #r2trim = "results/{sample}/misc/trimmed_reads/trimmed_reads2.fastq.gz",
        trim_finish = "results/{prefix}/{sample}/misc/trimmed_reads/perSVade_finished_file.txt",
    output:
        aligned_reads = "results/{prefix}/{sample}/misc/aligned_reads/aligned_reads.bam.sorted",
    params:
        r1trim = "results/{prefix}/{sample}/misc/trimmed_reads/trimmed_reads1.fastq.gz",
        r2trim = "results/{prefix}/{sample}/misc/trimmed_reads/trimmed_reads2.fastq.gz",
        outdir = "results/{prefix}/{sample}/misc/aligned_reads/",
        trim_outdir = "results/{prefix}/{sample}/misc/trimmed_reads/",
        ref = REFERENCE,
        fraction_available_mem = lambda wildcards, resources: round((resources.mem_mb/1000)/175 - 0.005,2),
    resources:
        mem_mb=10000,
        runtime=900,
    threads: 4
    priority: 100
    singularity:
        "docker://mikischikora/persvade:v1.02.6"
    shell:
        """
        set +u
        source /opt/conda/etc/profile.d/conda.sh
        conda activate perSVade_env
        python /perSVade/scripts/perSVade align_reads -f1 {params.r1trim} -f2 {params.r2trim} --ref {params.ref} -o {params.outdir} \
        --replace --verbose --fraction_available_mem {params.fraction_available_mem} --threads {threads}
        rm results/{wildcards.prefix}/{wildcards.sample}/misc/trimmed_reads/trimmed_reads*.zip
        rm results/{wildcards.prefix}/{wildcards.sample}/misc/trimmed_reads/trimmed_reads*.gz
        set -u
        """


# Step 3 will not be run normally!
# To manually activate this step, remove the commented input 'check_file' in step 4.
# If possible, this rule should be run with a single high-quality sample before running the full pipeline
# The output from a single run will be placed in the reference/parameter_optimization/ directory and can be re-used in future runs
# Note that this uses random SV simulations (as per the persvade FAQ) and assumes a haploid genome
rule step3_optimize_params:
    input:
        repeats = REF_DIR + 'repeat_inference/combined_repeats.tab',
        aligned_reads = "results/{prefix}/{sample}/misc/aligned_reads/aligned_reads.bam.sorted",
    output:
        check_file = "results/{prefix}/{sample}/misc/parameter_optimization/optp.txt",
    params:
        opt_params_dir = REF_DIR + 'parameter_optimization/',
        opt_params = REF_DIR + 'parameter_optimization/optimized_parameters.json',
        note_file = REF_DIR + 'parameter_optimization/{sample}_used_for_optimization.txt',
        ref = REFERENCE,
        mito_chr = config["mito_chr"],
        fraction_available_mem = lambda wildcards, resources: round((resources.mem_mb/1000)/175 - 0.005,2),
    resources:
        mem_mb=15000,
        runtime=480,
    threads: 4
    singularity:
        "docker://mikischikora/persvade:v1.02.6"
    shell:
        """
        set +u
        if [ -e {params.opt_params} ];
        then
            touch {output.check_file}
        else
            source /opt/conda/etc/profile.d/conda.sh
            conda activate perSVade_env
            mkdir -p {params.opt_params_dir}
            python /perSVade/scripts/perSVade optimize_parameters --ref {params.ref} -o {params.opt_params_dir} -sbam {input.aligned_reads} --repeats_file {input.repeats} \
            --mitochondrial_chromosome {params.mito_chr} --regions_SVsimulations random --simulation_ploidies haploid \
            --replace --verbose --fraction_available_mem {params.fraction_available_mem} --threads {threads}
            touch {params.note_file}
            touch {output.check_file}
        fi
        set -u
        """

# This step will automatically use optimized parameter files found in the reference folder!
# Step 3 must be run at least once before this step will work.
rule step4_call_svs:
    input:
        repeats = REF_DIR + 'repeat_inference/combined_repeats.tab',
        aligned_reads = "results/{prefix}/{sample}/misc/aligned_reads/aligned_reads.bam.sorted",
        ### Remove the line below to run parameter optimization
        #check_file = "results/{prefix}/{sample}/misc/parameter_optimization/optp.txt",
    output:
        call_svs = "results/{prefix}/{sample}/call_SVs/perSVade_finished_file.txt",
    params:
        opt_params = REF_DIR + 'parameter_optimization/optimized_parameters.json',
        outdir = "results/{prefix}/{sample}/call_SVs/",
        ref = REFERENCE,
        mito_chr = config["mito_chr"],
        fraction_available_mem = lambda wildcards, resources: round((resources.mem_mb/1000)/175 - 0.005,2),
    resources:
        mem_mb=10000,
        runtime=480,
    threads: 4
    singularity:
        "docker://mikischikora/persvade:v1.02.6"
    shell:
        """
        set +u
        source /opt/conda/etc/profile.d/conda.sh
        conda activate perSVade_env
        python /perSVade/scripts/perSVade call_SVs --ref {params.ref} -o {params.outdir} -sbam {input.aligned_reads} \
        --SVcalling_parameters {params.opt_params} --repeats_file {input.repeats} --mitochondrial_chromosome {params.mito_chr} \
        --replace --verbose --fraction_available_mem {params.fraction_available_mem} --threads {threads}
        set -u
        """


rule step5_vaf:
    input:
        call_svs = "results/{prefix}/{sample}/call_SVs/perSVade_finished_file.txt",
    output:
        vaf_file = "results/{prefix}/{sample}/call_SVs/vaf_finished.txt",
    params:
        outdir = "results/{prefix}/{sample}/call_SVs/",
    resources:
        mem_mb=3000,
        runtime=20,
    threads: 1
    shell:
        """
        python scripts/calculate_vaf.py --input {params.outdir}
        """

rule step6_cnv:
    input:
        aligned_reads = "results/{prefix}/{sample}/misc/aligned_reads/aligned_reads.bam.sorted",
    output:
        call_cnvs = "results/{prefix}/{sample}/call_CNVs/final_CNVcalling.tab",
    params:
        outdir = "results/{prefix}/{sample}/call_CNVs/",
        ref = REFERENCE,
        mito_chr = config["mito_chr"],
        fraction_available_mem = lambda wildcards, resources: round((resources.mem_mb/1000)/175 - 0.005,2),
    resources:
        mem_mb=10000,
        runtime=60,
    threads: 4
    singularity:
        "docker://mikischikora/persvade:v1.02.6"
    shell:
        """
        set +u
        source /opt/conda/etc/profile.d/conda.sh
        conda activate perSVade_env
        python /perSVade/scripts/perSVade call_CNVs --ref {params.ref} -o {params.outdir} -sbam {input.aligned_reads} \
        --mitochondrial_chromosome {params.mito_chr} --ploidy 1 --cnv_calling_algs HMMcopy,AneuFinder --window_size_CNVcalling 500 \
        --replace --verbose --fraction_available_mem {params.fraction_available_mem} --threads {threads}
        set -u
        """


rule step7_integrate:
    input:
        repeats = REF_DIR + 'repeat_inference/combined_repeats.tab',
        aligned_reads = "results/{prefix}/{sample}/misc/aligned_reads/aligned_reads.bam.sorted",
        vaf_file = "results/{prefix}/{sample}/call_SVs/vaf_finished.txt",
        call_cnvs = "results/{prefix}/{sample}/call_CNVs/final_CNVcalling.tab",
    output:
        combined_variants = "results/{prefix}/{sample}/integrated_SV_CNV_calls/SV_and_CNV_variant_calling.vcf",
    params:
        outdir = "results/{prefix}/{sample}/integrated_SV_CNV_calls/",
        cnv_dir = "results/{prefix}/{sample}/call_CNVs/",
        sv_dir = "results/{prefix}/{sample}/call_SVs/",
        ref = REFERENCE,
        mito_chr = config["mito_chr"],
        fraction_available_mem = lambda wildcards, resources: round((resources.mem_mb/1000)/175 - 0.005,2),
    resources:
        mem_mb=5000,
        runtime=30,
    threads: 1
    singularity:
        "docker://mikischikora/persvade:v1.02.6"
    shell:
        """
        set +u
        source /opt/conda/etc/profile.d/conda.sh
        conda activate perSVade_env
        python /perSVade/scripts/perSVade integrate_SV_CNV_calls -o {params.outdir} --ref {params.ref} \
        --mitochondrial_chromosome {params.mito_chr} --ploidy 1 -sbam {input.aligned_reads} \
        --outdir_callSVs {params.sv_dir} --outdir_callCNVs {params.cnv_dir} --repeats_file {input.repeats} \
        --replace --verbose --fraction_available_mem {params.fraction_available_mem} --threads {threads}
        set -u
        """

# this requires mitochondrial and genomic codon tables
# these are currently set to 4 and 12 respectively for c auris
rule step8_annotate:
    input:
        combined_variants = "results/{prefix}/{sample}/integrated_SV_CNV_calls/SV_and_CNV_variant_calling.vcf",
    output:
        annotated_variants = "results/{prefix}/{sample}/annotate_SVs/annotated_variants.tab",
    params:
        outdir = "results/{prefix}/{sample}/annotate_SVs/",
        ref = REFERENCE,
        ref_gff = REF_GFF,
        mito_chr = config["mito_chr"],
        fraction_available_mem = lambda wildcards, resources: round((resources.mem_mb/1000)/175 - 0.005,2),
    resources:
        mem_mb=5000,
        runtime=30,
    threads: 1
    singularity:
        "docker://mikischikora/persvade:v1.02.6"
    shell:
        """
        set +u
        source /opt/conda/etc/profile.d/conda.sh
        conda activate perSVade_env
        python /perSVade/scripts/perSVade annotate_SVs -o {params.outdir} --ref {params.ref} \
        --mitochondrial_chromosome {params.mito_chr} -gff {params.ref_gff} -mcode 4 -gcode 12 \
        --SV_CNV_vcf {input.combined_variants} \
        --replace --verbose --fraction_available_mem {params.fraction_available_mem} --threads {threads}
        set -u
        """

