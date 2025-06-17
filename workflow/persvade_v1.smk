import subprocess
import pandas as pd
import os

configfile: "config/config.yaml"

samples_df = pd.read_csv(config["samples"])
SAMPLE = list(samples_df['sample_id'])

# prepare reference directory
# note that difference reference genomes are not supported yet - the reference must be the same for all samples!
def prep_reference(ref_file):
    ref_name = ref_file.split('/')[-1].split('.fa')[0]
    ref_path = f'results/reference/{ref_name}/'
    new_ref_file = f'results/reference/{ref_name}/{ref_name}.fasta'
    print(f'The reference used for ALL samples in this run is {ref_name}, located at {ref_file}')
    if not os.path.exists(ref_path):
        subprocess.call(['mkdir','-p',ref_path])
        subprocess.call(['cp',ref_file,new_ref_file])
        print(f'A new reference directory has been created at {ref_path}')
    else:
        print('Using existing reference files at this location!')
    print('Remember to check the reference genome mitochondrial chromosome location in config/config.yaml')
    return([new_ref_file,ref_path])

# this will copy the provided reference to a new directory (assuming it's not already there)
REFERENCE,REF_DIR = prep_reference(config["reference"])

rule all:
    input:
        call_svs = expand("results/{sample}/call_SVs/perSVade_finished_vaf_file.txt",sample=SAMPLE),

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


rule step1_qc:
    input:
        r1 = config["short_reads"] + "/" + "{sample}_R1.fastq.gz",
        r2 = config["short_reads"] + "/" + "{sample}_R2.fastq.gz",
    output:
        r1trim = "results/{sample}/misc/trimmed_reads/trimmed_reads1.fastq.gz",
        r2trim = "results/{sample}/misc/trimmed_reads/trimmed_reads2.fastq.gz",
    params:
        outdir = "results/{sample}/misc/trimmed_reads/",
        fraction_available_mem = lambda wildcards, resources: round((resources.mem_mb/1000)/175 - 0.005,2),
    resources:
        mem_mb=5000,
        runtime=40,
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
        r1trim = "results/{sample}/misc/trimmed_reads/trimmed_reads1.fastq.gz",
        r2trim = "results/{sample}/misc/trimmed_reads/trimmed_reads2.fastq.gz",
    output:
        aligned_reads = "results/{sample}/misc/aligned_reads/aligned_reads.bam.sorted",
    params:
        outdir = "results/{sample}/misc/aligned_reads/",
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

# If possible, this rule should be run with a single high-quality sample before running the full pipeline
# The output from a single run will be placed in the reference/parameter_optimization/ directory and can be re-used in future runs
# Note that this uses random SV simulations (as per the persvade FAQ) and assumes a haploid genome
rule step3_optimize_params:
    input:
        repeats = REF_DIR + 'repeat_inference/combined_repeats.tab',
        aligned_reads = "results/{sample}/misc/aligned_reads/aligned_reads.bam.sorted",
    output:
        check_file = "results/{sample}/misc/parameter_optimization/optp.txt",
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
rule step4_call_svs:
    input:
        repeats = REF_DIR + 'repeat_inference/combined_repeats.tab',
        aligned_reads = "results/{sample}/misc/aligned_reads/aligned_reads.bam.sorted",
        check_file = "results/{sample}/misc/parameter_optimization/optp.txt",
    output:
        call_svs = "results/{sample}/call_SVs/perSVade_finished_file.txt",
    params:
        opt_params = REF_DIR + 'parameter_optimization/optimized_parameters.json',
        outdir = "results/{sample}/call_SVs/",
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


rule step5_vaf:
    input:
        call_svs = "results/{sample}/call_SVs/perSVade_finished_file.txt",
    output:
        vaf_file = "results/{sample}/call_SVs/perSVade_finished_vaf_file.txt",
    params:
        outdir = "results/{sample}/call_SVs/",
        fraction_available_mem = lambda wildcards, resources: round((resources.mem_mb/1000)/175 - 0.005,2),
    resources:
        mem_mb=3000,
        runtime=20,
    threads: 1
    shell:
        """
        python /scripts/calculate_vaf.py --input {params.outdir}
        """


