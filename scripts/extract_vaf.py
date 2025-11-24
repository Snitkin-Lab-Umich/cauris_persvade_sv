import os
import argparse
import pandas as pd
import gffutils as gff

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
            scaffold1,start1,end1 = line_data[0],line_data[1],line_data[2]
            if sv_type in ['insertions','translocations']:
                scaffold2,start2,end2 = line_data[3],line_data[4],line_data[5]
            else:
                scaffold2,start2,end2 = 'NA','NA','NA'
            vaf_data = [float(x) for x in vaf_data.split('+')]
            d.append([isolatename,sv_type,str(round(min(vaf_data),2)),str(round(max(vaf_data),2)),scaffold1,start1,end1,scaffold2,start2,end2])
    final_df = pd.DataFrame(d,columns = ['isolate','sv_type','min_vaf','max_vaf','scaffold1','start1','end1','scaffold2','start2','end2'])
    return(final_df)

def summarize_vaf(dirpath):
    """return a single dataframe with all VAF values when given the expected output directory"""
    sv_columns = ['isolate','sv_type','min_vaf','max_vaf','scaffold1','start1','end1','scaffold2','start2','end2']
    # make a single df to store data from all isolates
    outdata = pd.DataFrame([],columns = sv_columns)
    for isolatename in os.listdir(dirpath):
        # make a small df to store data for a single isolate
        isolate_outdata = pd.DataFrame([],columns = sv_columns)
        sv_path = os.path.join(dirpath, isolatename, 'call_SVs')
        if not os.path.isdir(sv_path):
            print(f'Could not locate SV calling directory for isolate {isolatename}')
            continue
        for f in os.listdir(sv_path):
            if f.endswith('_vaf.tab'):
                filepath = os.path.join(sv_path, f)
                vaf_df = get_vaf_values(filepath, isolatename)
                isolate_outdata = pd.concat([isolate_outdata, vaf_df], ignore_index=True)
        # save the isolate data as its own csv
        isolate_outdata.to_csv(os.path.join(sv_path, 'vaf_summary.csv'), index=False)
        # add the isolate data to the overall dataframe
        outdata = pd.concat([outdata, isolate_outdata], ignore_index=True)
    return(outdata)

def add_gene_data(vaf_df, gff_path):
    # take a dataframe containing genome coordinates
    # add a column containing all cds features from the gff that overlap with those coordinates
    db = gff.create_db(gff_path, dbfn=':memory:', force=True, keep_order=True, merge_strategy='merge', sort_attribute_values=True)
    vaf_df['genes_overlapping_sv'] = ''
    for r in range(len(vaf_df)):
        sv_type = vaf_df.loc[r,'sv_type']
        scaffold1 = vaf_df.loc[r,'scaffold1']
        start1 = int(float(vaf_df.loc[r,'start1']))
        end1 = int(float(vaf_df.loc[r,'end1']))
        overlapping_genes = set()
        for feature in db.region(seqid=scaffold1, start=start1, end=end1, completely_within=False):
            if feature.featuretype == 'CDS':
                overlapping_genes.add(feature.attributes.get('locus_tag',[feature.id])[0])
        if sv_type in ['insertions','translocations']:
            scaffold2 = vaf_df.loc[r,'scaffold2']
            start2 = int(float(vaf_df.loc[r,'start2']))
            end2 = int(float(vaf_df.loc[r,'end2']))
            for feature in db.region(seqid=scaffold2, start=start2, end=end2, completely_within=False):
                if feature.featuretype == 'CDS':
                    overlapping_genes.add(feature.attributes.get('locus_tag',[feature.id])[0])
        vaf_df.at[r,'genes_overlapping_sv'] = ';'.join(sorted(overlapping_genes)) if overlapping_genes else 'NA'
    return(vaf_df)



def main():
    # define all args
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--input','-i',type=str,
        help='''Provide a path to a persvade snakemake output directory.''',
        required=True
        )
    parser.add_argument(
        '--gff','-g',type=str,
        help='''Provide a path to the reference gff used for this persvade run.''',
        required=True
        )
    parser.add_argument(
        '--output','-o',type=str,
        help='''Provide a file path for the output csv.''',
        default='vaf_summary.csv'
        )
    args = parser.parse_args()
    outdata = summarize_vaf(args.input)
    outdata2 = add_gene_data(outdata, args.gff)
    outdata2.to_csv(args.output, index=False)


if __name__ == "__main__":
    main()



