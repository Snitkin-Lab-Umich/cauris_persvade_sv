import os
import argparse
import subprocess

def read_gridss_info(vcf):
    """Read the INFO field from a GRIDSS VCF file and return a dictionary."""
    gridss_dict = {}
    with open(vcf, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            gridss_name = parts[2]
            gridss_info = parts[7]
            for item in gridss_info.split(';'):
                if '=' in item:
                    key, value = item.split('=', 1)
                    if gridss_name not in gridss_dict:
                        gridss_dict[gridss_name] = {key: value}
                    else:
                        gridss_dict[gridss_name][key] = value
    return gridss_dict

def add_vaf_to_sv_file(sv_file, gridss_info, gridss_threshold_single=0.8, gridss_threshold_all=0.2, fragsize=1000):
    """Read a df of structural variants from persvade and write a new version that has a VAF column."""
    # the default VAF thesholds are from https://www.nature.com/articles/s41564-023-01547-z
    sv_type = sv_file.split('/')[-1].split('.')[0]
    new_sv_file = sv_file.replace('.tab', '_vaf.tab')
    with open(sv_file, 'r') as fh_in, open(new_sv_file, 'w') as fh_out:
        # add twp extra columns to the header for the VAF values and the threshold
        _ = fh_out.write(fh_in.readline().replace('\n', '\tVAF\tQC_EVAL\n'))
        for line in fh_in:
            line = line.strip().split('\t')
            # extract the GRIDSS names and SV length, depending on the SV type
            if sv_type in ['deletions','inversions','tandemDuplications']:
                # these have TWO columns for start and end
                sv_length = abs(int(line[1]) - int(line[2]))
                gridss_list = line[3].split('+')
            elif sv_type in ['insertions','translocations']:
                # these have FOUR columns for start and end
                # both lengths should be comparable, so the shorter one is taken
                sv_length1 = abs(int(line[1]) - int(line[2]))
                sv_length2 = abs(int(line[4]) - int(line[5]))
                if sv_length1/sv_length2 > 1.5 or sv_length2/sv_length1 > 1.5:
                    # if the lengths differ too much, print a warning
                    print(f'Warning: Differing lengths for {sv_type} found: {sv_length1} vs {sv_length2}. Using the shorter one.')
                sv_length = min(sv_length1, sv_length2)
                gridss_list = line[7].split('+')
            # calculate the VAF for each GRIDSS name in the list
            vaf_scores = []
            for gridss_name in gridss_list:
                if gridss_name in gridss_info:
                    vf = int(gridss_info[gridss_name]['VF'])
                    ref = int(gridss_info[gridss_name]['REF'])
                    refpair = int(gridss_info[gridss_name]['REFPAIR'])
                    if sv_length >= fragsize:
                        vaf = round(vf / (vf + ref + refpair),3)
                    else:
                        vaf = round(vf / (vf + ref),3)
                    vaf_scores.append(vaf)
                else:
                    print(f'Warning: GRIDSS name {gridss_name} not found in GRIDSS info. Setting VAF to 0.')
                    vaf_scores.append(0.0)
            # determine if the VAF scores are above the threshold
            check1 = any([v >= gridss_threshold_single for v in vaf_scores])
            check2 = all([v >= gridss_threshold_all for v in vaf_scores])
            if check1 and check2:
                qc_eval = 'PASS'
            else:
                qc_eval = 'FAIL'
            # write the line to the new file
            vaf_scores_str = [str(v) for v in vaf_scores]
            fh_out.write('\t'.join(line) + f'\t{"+".join(vaf_scores_str)}\t{qc_eval}\n')

# remember to read the vcf info only once, not for each SV
def sv_file_iter(sv_file_list, vcf_file):
    # start by reading in the gridss info
    gridss_info = read_gridss_info(vcf_file)
    # for each SV file, read the SVs and calculate VAF
    for sv_file in sv_file_list:
        sv_file_path = os.path.abspath(sv_file)
        add_vaf_to_sv_file(sv_file_path, gridss_info)
        print(f'Added VAF to {sv_file}.')

def main():
    # define all args
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--input','-i',type=str,
        help='''Provide a path to a persvade call_SVs output directory. These should all end in .tab.''',
        required=True
        )
    parser.add_argument(
        '--gridss','-g',type=str,
        help='''Provide an alternate path to a GRIDSS VCF file. Without this option, the script 
        will look for gridss_output.filt.vcf in the persvade call_SVs directory provided above.''',
        default=None
        )
    args = parser.parse_args()
    if not args.input.endswith('/'):
        args.input+='/'
    sv_file_list = []
    for f in os.listdir(args.input):
        if f.endswith('.tab') and any([x in f for x in ['deletions','inversions','tandemDuplications','insertions','translocations']]):
            if '_vaf.tab' not in f:
                # avoid reprocessing files
                sv_file_list.append(os.path.join(args.input, f))
    if args.gridss is None:
        args.gridss = os.path.join(args.input, 'gridss_output.filt.vcf')
    if not os.path.exists(args.gridss):
        raise FileNotFoundError(f'GRIDSS VCF file {args.gridss} not found. Please provide a valid path.')
    sv_file_iter(sv_file_list,args.gridss)
    subprocess.call(['touch', os.path.join(args.input, 'perSVade_finished_vaf_file.txt')])

if __name__ == "__main__":
    main()



