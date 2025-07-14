import os
import argparse
import pandas as pd

def get_vaf_values(filepath,isolatename):
    """Return the highest and lowest VAF scores from a given VAF SV file"""
    sv_type = filepath.split('/')[-1].split('_vaf.tab')[0]
    if sv_type in ['deletions','inversions','tandemDuplications']:
        vaf_position = 4
    if sv_type in ['insertions','translocations']:
        vaf_position = 8
    d = []
    with open(filepath,'r') as fh:
        _ = next(fh)
        for line in fh:
            line_data = line.strip().split('\t')
            vaf_data = line_data[vaf_position]
            scaffold = line_data[0]
            vaf_data = [float(x) for x in vaf_data.split('+')]
            d.append([isolatename,sv_type,str(round(min(vaf_data),2)),str(round(max(vaf_data),2)),scaffold])
    return(pd.DataFrame(d,columns = ['isolate','sv_type','min_vaf','max_vaf','scaffold']))

def summarize_vaf(dirpath):
    """return a single dataframe with all VAF values when given the expected output directory"""
    outdata = pd.DataFrame([],columns = ['isolate','sv_type','min_vaf','max_vaf','scaffold'])
    for isolatename in os.listdir(dirpath):
        sv_path = os.path.join(dirpath, isolatename, 'call_SVs')
        if not os.path.isdir(sv_path):
            continue
        for f in os.listdir(sv_path):
            if f.endswith('_vaf.tab'):
                filepath = os.path.join(sv_path, f)
                vaf_df = get_vaf_values(filepath, isolatename)
                outdata = pd.concat([outdata, vaf_df], ignore_index=True)
    return(outdata)

def main():
    # define all args
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--input','-i',type=str,
        help='''Provide a path to a persvade snakemake output directory.''',
        required=True
        )
    parser.add_argument(
        '--output','-o',type=str,
        help='''Provide a file path for the output csv.''',
        default='vaf_summary.csv'
        )
    args = parser.parse_args()
    outdata = summarize_vaf(args.input)
    outdata.to_csv(args.output, index=False)


if __name__ == "__main__":
    main()



