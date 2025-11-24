import os
import argparse
import subprocess
import pandas as pd

def search_samples(input_csv, trimmed_reads_dirlist, raw_reads_dirlist, out_dir):
    # search through an input csv containing sample names
    # for each sample, first look in the trimmed reads directories to see if the expected file is present
    # if it's missing, look in the raw read directories and add it to the output list
    # finally, write a new csv with all samples that only have raw reads
    df = pd.read_csv(input_csv,sep='\t')
    samples_list = df['Sample'].tolist()
    # make all paths
    out_csv_path = os.path.join(out_dir, 'samples_to_trim.csv')
    trim_dir_out = os.path.join(out_dir, 'trimmed_reads')
    raw_dir_out = os.path.join(out_dir, 'raw_reads')
    # make all output directories
    os.makedirs(trim_dir_out, exist_ok=True)
    os.makedirs(raw_dir_out, exist_ok=True)
    # iterate through all samples
    to_trim_list = []
    for sampl in samples_list:
        found_trimmed,found_raw = False,False
        # first, search all trimmed reads directories
        r1_trim_file = f'{sampl}_R1_trim_paired.fastq.gz'
        r2_trim_file = f'{sampl}_R2_trim_paired.fastq.gz'
        for trim_dir in trimmed_reads_dirlist:
            trim_r1_path = os.path.join(trim_dir, r1_trim_file)
            trim_r2_path = os.path.join(trim_dir, r2_trim_file)
            if os.path.isfile(trim_r1_path) and os.path.isfile(trim_r2_path):
                # copy the trimmed reads to the output directory
                subprocess.run(['cp', trim_r1_path, os.path.join(trim_dir_out, r1_trim_file)])
                subprocess.run(['cp', trim_r2_path, os.path.join(trim_dir_out, r2_trim_file)])
                found_trimmed = True
                break
        # if none are found, search raw reads directories
        if not found_trimmed:
            r1_raw_file = f'{sampl}_R1.fastq.gz'
            r2_raw_file = f'{sampl}_R2.fastq.gz'
            for raw_dir in raw_reads_dirlist:
                raw_r1_path = os.path.join(raw_dir, r1_raw_file)
                raw_r2_path = os.path.join(raw_dir, r2_raw_file)
                if os.path.isfile(raw_r1_path) and os.path.isfile(raw_r2_path):
                    # copy the raw reads to the output directory
                    subprocess.run(['cp', raw_r1_path, os.path.join(raw_dir_out, r1_raw_file)])
                    subprocess.run(['cp', raw_r2_path, os.path.join(raw_dir_out, r2_raw_file)])
                    found_raw = True
                    to_trim_list.append(sampl)
                    break
        # all samples should have either trimmed or raw reads
        if not found_trimmed and not found_raw:
            print(f'Error: Sample {sampl} not found in any provided directories.')
            quit(1)
    # write the new csv file
    with open(out_csv_path, 'w') as fh_out:
        fh_out.write('Sample\n')
        for sampl in to_trim_list:
            fh_out.write(f'{sampl}\n')


def main():
    # define all args
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--input','-i',type=str,
        help='''Provide an input samples.csv file.''',
        required=True
        )
    parser.add_argument(
        '--trimmed_reads_dirs','-tdir',type=str,nargs='+',
        help='''Provide a list of paths to directories containing trimmed reads.''',
        default=None
        )
    parser.add_argument(
        '--raw_reads_dirs','-rdir',type=str,nargs='+',
        help='''Provide a list of paths to directories containing raw reads.''',
        default=None
        )
    parser.add_argument(
        '--output','-o',type=str,
        help='''Provide a path to the output directory.''',
        default=None
        )
    args = parser.parse_args()
    search_samples(
        input_csv=args.input,
        trimmed_reads_dirlist=args.trimmed_reads_dirs,
        raw_reads_dirlist=args.raw_reads_dirs,
        out_dir=args.output
    )

if __name__ == "__main__":
    main()



