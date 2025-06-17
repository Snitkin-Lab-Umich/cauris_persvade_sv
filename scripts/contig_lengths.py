
import argparse


def count_contig_len(input_file):
    with open(input_file,'r') as fh:
        t = None
        data_list = []
        for line in fh:
            line = line.strip()
            if line[0] == '>':
                if t is not None:
                    data_list.append(t)
                contig_name = line.split('>')[1]
                # remove anything after a space to keep contig names consistent with mummer and BLAST
                contig_name = contig_name.split(' ')[0]
                t = [contig_name,0]
            else:
                t[1]+=len(line)
    data_list.append(t)
    for el in data_list:
        print(f'{el[0]}: {el[1]}')



def main():
    # define all args
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--input','-i',type=str,
        help='''Provide a path to an assembly (.fasta).''',
        required=True
        )
    args = parser.parse_args()
    count_contig_len(args.input)

if __name__ == "__main__":
    main()



