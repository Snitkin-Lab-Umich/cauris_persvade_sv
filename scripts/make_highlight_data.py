import os
import argparse
import pandas as pd


def make_hdata(input_path, reference_path):
    hdata = parse_single(input_path + 'deletions.tab')
    hdata = pd.concat([hdata,parse_single(input_path + 'inversions.tab')])
    hdata = pd.concat([hdata,parse_single(input_path + 'tandemDuplications.tab')])
    hdata = pd.concat([hdata,parse_double(input_path + 'insertions.tab')])
    #hdata = pd.concat([hdata,parse_double(input_path + 'translocations.tab')])
    hdata['type'] = ['subject'] * len(hdata)
    refname = reference_path.split('/')[-1].split('.')[0]
    hdata['name'] = [refname] * len(hdata)
    hdata = hdata[['type','name','contig','start','end']]
    hdata.to_csv('highlight_data.tsv',sep='\t',index=False)


def parse_single(infile):
    del_data = pd.read_csv(infile,sep='\t')
    del_data_highlight = pd.DataFrame({'contig': del_data['Chr'],'start': del_data['Start'],'end': del_data['End']})
    return(del_data_highlight)

def parse_double(infile):
    del_data = pd.read_csv(infile,sep='\t')
    del_data_highlight_ChrA = pd.DataFrame({'contig': del_data['ChrA'],'start': del_data['StartA'],'end': del_data['EndA']})
    del_data_highlight_ChrB = pd.DataFrame({'contig': del_data['ChrB'],'start': del_data['StartB'],'end': del_data['EndB']})
    return(pd.concat([del_data_highlight_ChrA,del_data_highlight_ChrB]))



def main():
    # define all args
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--input','-i',type=str,
        help='''Provide a path to a directory with the output of perSVade callSVs. This should contain deletions.tab, insertions.tab, 
        inversions.tab, translocations.tab, and tandemDuplications.tab.''',
        required=True
        )
    parser.add_argument(
        '--reference','-r',type=str,
        help='''Provide a path to the reference file used when running perSVade. Use this same file with the dotplot script.''',
        required=True
        )
    args = parser.parse_args()
    if not args.input.endswith('/'):
        args.input+='/'
    for p in [args.input,args.reference]:
        if not os.path.exists(p):
            print(f'Could not locate {p}')
            quit(1)
    make_hdata(args.input,args.reference)


if __name__ == "__main__":
    main()



