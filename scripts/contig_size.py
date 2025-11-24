from Bio import SeqIO
import argparse
import os

def contig_lengths(fasta_path):
    with open(fname,'r') as fh:
    for record in SeqIO.parse(fh,'fasta'):
        print(f'{record.id}: {str(len(record.seq))}')


def main():
    # define all args
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--input','-i',type=str,
        help='''Provide a path to a fasta.''',
        required=True
        )
    args = parser.parse_args()
    if not os.path.isfile(args.input):
        print(f'Could not locate {args.input}')
        quit(1)
    contig_lengths(args.input)

if __name__ == "__main__":
    main()



