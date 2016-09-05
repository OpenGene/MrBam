#!/usr/bin/env python3

import argparse
import bam
import aggregate

parser = argparse.ArgumentParser()

parser.add_argument('query', help="vcf file contains mutations to query")
parser.add_argument('reads', help="bam file contains origin reads info. There must be a corresponding .bai file in the same directory")
parser.add_argument('--quality', '-q', type=int, default=20, help="drop bases whose qulity is less than this")
parser.add_argument('--ref', '-r', help="reference file. If not provided, consider the mode of reads as the ref")
parser.add_argument('--verbos', '-v', help="output debug info")

o = parser.parse_args()

with open(o.query + ".readsinfo.txt", "w") as fout:
    fout.write('#Chr' + '\t' + 'Position' + '\t' + 'Unique_Pairs_Support_Alt' + '\t' + 'Unique_Single_Support_Alt' + '\t')
    fout.write('N_Reads' + '\t' + 'N_Error' + '\t' + 'N_Low_Quality' + '\t' + 'N_Inconsistence_Pairs' + '\n')

    with open(o.query) as fin:
        for line in fin:
            if line.startswith('#'):
                continue

            p = ':'.join(line.split('\t')[:2])

            unique_pairs, unique_single, nsum, nerror, nlowq, ninconsis = aggregate.aggregate_reads(o, bam.get_reads(o, p))

            d = {'A': 0, 'T': 0, 'C': 0, 'G': 0}

            for x in (unique_single, unique_pairs):
                for startpos, endpos, base in x:
                    d[base] += 1

            refbase = max(d, key=lambda x: d[x])

            unique_pairs_support_alt  = sum(1 for x in unique_pairs if x[2]!=refbase)
            unique_single_support_alt = sum(1 for x in unique_single if x[2]!=refbase)

            fout.write(p.replace(':', '\t') + '\t' + str(unique_pairs_support_alt) + '\t' + str(unique_single_support_alt) + '\t')
            fout.write(str(nsum) + '\t' + str(nerror) + '\t' + str(nlowq) + '\t' + str(ninconsis) + '\n')
