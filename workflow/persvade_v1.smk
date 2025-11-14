import subprocess
import os
import pandas as pd

configfile: "config/config.yaml"

samples_df = pd.read_csv(config['samples'],sep = '\t')
SAMPLE = samples_df['Sample'].tolist()

PREFIX = config["prefix"]

# prepare reference directory
def prep_reference(ref_fasta, ref_gff):
    ref_name = ref_fasta.split('/')[-1].split('.fa')[0]
    ref_path = f'results/reference/{ref_name}/'
    new_ref_fasta = f'results/reference/{ref_name}/{ref_name}.fasta'
    new_ref_gff = f'results/reference/{ref_name}/{ref_name}.gff'
    print(f'The reference used for this run is {ref_name}, located at {ref_fasta}')
    if not os.path.exists(new_ref_path) or not os.path.exists(new_ref_gff):
        subprocess.call(['mkdir','-p',ref_path])
        subprocess.call(['cp',ref_fasta,new_ref_fasta])
        subprocess.call(['cp',ref_gff,new_ref_gff])
        print(f'A new reference directory has been created at {ref_path}')
    else:
        print('Using existing reference files at this location!')
    print('Remember to check the mitochondrial chromosome location in config/config.yaml')
    return([new_ref_fasta,new_ref_gff,ref_path])

# this will copy the provided reference to a new directory (assuming it's not already there)
REFERENCE,REF_GFF,REF_DIR = prep_reference(config["reference_fasta"],config["reference_gff"])

rule all:
    input:
        annotate = expand("results/{prefix}/{sample}/annotate_SVs/annotated_variants.tab",sample=SAMPLE),

# the qc step should not be necessary, since we run trimmomatic and qc during the assembly and annotation
# rule step1_qc:
#     input:
#         r1 = config["short_reads"] + "/" + "{sample}_R1_trim_paired.fastq.gz",
#         r2 = config["short_reads"] + "/" + "{sample}_R2_trim_paired.fastq.gz",
#     output:
#         r1trim = "results/{prefix}/{sample}/misc/trimmed_reads/trimmed_reads1.fastq.gz",
#         r2trim = "results/{prefix}/{sample}/misc/trimmed_reads/trimmed_reads2.fastq.gz",
#     params:
#         outdir = "results/{prefix}/{sample}/misc/trimmed_reads/",
#         fraction_available_mem = lambda wildcards, resources: round((resources.mem_mb/1000)/175 - 0.005,2),
#     resources:
#         mem_mb=5000,
#         runtime=40,
#     threads: 1
#     singularity:
#         "docker://mikischikora/persvade:v1.02.6"
#     shell:
#         """
#         set +u
#         source /opt/conda/etc/profile.d/conda.sh
#         conda activate perSVade_env
#         python /perSVade/scripts/perSVade trim_reads_and_QC -f1 {input.r1} -f2 {input.r2} -o {params.outdir} \
#         --replace --verbose --fraction_available_mem {params.fraction_available_mem} --threads {threads}
#         set -u
#         """

rule step2_align:
    input:
        r1 = config["trimmed_reads"] + "/" + "{sample}_R1_trim_paired.fastq.gz",
        r2 = config["trimmed_reads"] + "/" + "{sample}_R2_trim_paired.fastq.gz",
    output:
        aligned_reads = "results/{prefix}/{sample}/misc/aligned_reads/aligned_reads.bam.sorted",
    params:
        outdir = "results/{prefix}/{sample}/misc/aligned_reads/",
        ref = REFERENCE,
        fraction_available_mem = lambda wildcards, resources: round((resources.mem_mb/1000)/175 - 0.005,2),
    resources:
        mem_mb=15000,
        runtime=360,
    threads: 4
    singularity:
        "docker://mikischikora/persvade:v1.02.6"
    shell:
        """
        set +u
        source /opt/conda/etc/profile.d/conda.sh
        conda activate perSVade_env
        python /perSVade/scripts/perSVade align_reads -f1 {input.r1trim} -f2 {input.r2trim} --ref {params.ref} -o {params.outdir} \
        --replace --verbose --fraction_available_mem {params.fraction_available_mem} --threads {threads}
        set -u
        """

# Repeat inference only needs to be run once per reference genome
rule step3_repeats:
    input:
        ref = REFERENCE,
    output:
        repeats = REF_DIR + 'repeat_inference/combined_repeats.tab',
    params:
        outdir = REF_DIR + 'repeat_inference/',
        fraction_available_mem = lambda wildcards, resources: round((resources.mem_mb/1000)/175 - 0.005,2),
    resources:
        mem_mb=15000,
        runtime=360,
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


# If possible, this rule should be run with a single high-quality sample before running the full pipeline
# The output from a single run will be placed in the reference/parameter_optimization/ directory and can be re-used in future runs
# Note that this uses random SV simulations (as per the persvade FAQ) and assumes a haploid genome
rule step4_optimize_params:
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
rule step5_1_call_svs:
    input:
        repeats = REF_DIR + 'repeat_inference/combined_repeats.tab',
        aligned_reads = "results/{prefix}/{sample}/misc/aligned_reads/aligned_reads.bam.sorted",
        check_file = "results/{prefix}/{sample}/misc/parameter_optimization/optp.txt",
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
        runtime=60,
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


rule step5_2_vaf:
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


# the cnv calling step does not require the repeat file for some reason
rule step6_call_cnvs:
    input:
        repeats = REF_DIR + 'repeat_inference/combined_repeats.tab',
        aligned_reads = "results/{prefix}/{sample}/misc/aligned_reads/aligned_reads.bam.sorted",
        check_file = "results/{prefix}/{sample}/misc/parameter_optimization/optp.txt",
    output:
        call_cnvs = "results/{prefix}/{sample}/call_CNVs/perSVade_finished_file.txt",
    params:
        outdir = "results/{prefix}/{sample}/call_CNVs/",
        ref = REFERENCE,
        mito_chr = config["mito_chr"],
        ploidy = config["ploidy"],
        cnv_algs = config["cnv_calling_algorithms"],
        window_size = config["cnv_calling_window_size"],
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
        --mitochondrial_chromosome {params.mito_chr} --ploidy {params.ploidy} \
        --cnv_calling_algs {params.cnv_algs} --window_size_CNVcalling {params.window_size} \
        --replace --verbose --fraction_available_mem {params.fraction_available_mem} --threads {threads}
        set -u
        """

rule step7_integrate:
    input:
        repeats = REF_DIR + 'repeat_inference/combined_repeats.tab',
        aligned_reads = "results/{prefix}/{sample}/misc/aligned_reads/aligned_reads.bam.sorted",
        check_file = "results/{prefix}/{sample}/misc/parameter_optimization/optp.txt",
        call_svs = "results/{prefix}/{sample}/call_SVs/perSVade_finished_file.txt",
        call_cnvs = "results/{prefix}/{sample}/call_CNVs/perSVade_finished_file.txt",
    output:
        integrate = "results/{prefix}/{sample}/integrated_SV_CNV_calls/SV_and_CNV_variant_calling.vcf",
    params:
        outdir = "results/{prefix}/{sample}/integrated_SV_CNV_calls/",
        ref = REFERENCE,
        mito_chr = config["mito_chr"],
        ploidy = config["ploidy"],
        sv_dir = "results/{prefix}/{sample}/call_SVs/",
        cnv_dir = "results/{prefix}/{sample}/call_CNVs/",
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
        python /perSVade/scripts/perSVade integrate_SV_CNV_calls --ref {params.ref} -o {params.outdir} -sbam {input.aligned_reads} \
        --outdir_callSVs {params.sv_dir} --outdir_callCNVs {params.cnv_dir} --repeats_file {input.repeats} \
        --mitochondrial_chromosome {params.mito_chr} --ploidy {params.ploidy} \
        --replace --verbose --fraction_available_mem {params.fraction_available_mem} --threads {threads}
        set -u
        """

rule step8_annotate:
    input:
        repeats = REF_DIR + 'repeat_inference/combined_repeats.tab',
        aligned_reads = "results/{prefix}/{sample}/misc/aligned_reads/aligned_reads.bam.sorted",
        check_file = "results/{prefix}/{sample}/misc/parameter_optimization/optp.txt",
        integrate = "results/{prefix}/{sample}/integrated_SV_CNV_calls/SV_and_CNV_variant_calling.vcf",
    output:
        annotate = "results/{prefix}/{sample}/annotate_SVs/annotated_variants.tab",
    params:
        outdir = "results/{prefix}/{sample}/annotate_SVs/",
        ref = REFERENCE,
        ref_gff = REF_GFF,
        mito_chr = config["mito_chr"],
        gcode = config["gcode"],
        mcode = config["mcode"],
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
        python /perSVade/scripts/perSVade annotate_SVs --ref {params.ref} -gff {params.ref_gff} -o {params.outdir} \
        --SV_CNV_vcf {input.integrate} --gcode {params.gcode} --mcode {params.mcode} \
        --mitochondrial_chromosome {params.mito_chr} \
        --replace --verbose --fraction_available_mem {params.fraction_available_mem} --threads {threads}
        set -u
        """

