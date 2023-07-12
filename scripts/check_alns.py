#!/usr/bin/env python3

import argparse
import glob

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", type=str, required=True)
parser.add_argument("-o", "--output", type=str)


args = parser.parse_args()

aln_dir = args.input

def is_weird(m):
    with open(m) as infile:
        first_letter = (''.join(line[0] for line in infile if line))
        if first_letter == 'b':
            return True
        else:
            return False

if __name__ == "__main__":
    alns = glob.glob(aln_dir+'/*')
    for aln in alns:
        with open(aln, 'r+') as file:
            if is_weird(aln):
                lines = file.readlines()
                file.seek(0)
                for line in lines:
                    newline = line.replace('b', '').replace('\\n', '\n').replace('\'','')
                    file.write(newline)
                file.truncate()
            else:
                exit
    # with open(gene_file) as input:
    #     genes_dict = [read_treeline(line) for line in input]
    # print('Done loading genes!')
