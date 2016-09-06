#!/usr/bin/env python3

import argparse
import bam
import aggregate

parser = argparse.ArgumentParser()

parser.add_argument('query', help="vcf file contains mutations to query")
parser.add_argument('reads', help="bam file contains origin reads info. There must be a corresponding .bai file in the same directory")
parser.add_argument('--quality', '-q', type=int, default=20, help="drop bases whose qulity is less than this")
parser.add_argument('--ref', '-r', help="reference file. If not provided, consider the mode of reads as the ref")
parser.add_argument('--verbos', '-v', action='store_true', help="output debug info")

o = parser.parse_args()

with open(o.query + ".readsinfo.txt", "w") as fout:
    println = lambda *args: fout.write('\t'.join(map(str, args)) + '\n')

    println(
        '#Chr', 'Position',
        'Unique_Overlapped_Pairs_Support_Alt',
        'Unique_Non_Overlapped_Pairs_Support_Alt',
        'Unique_Single_Support_Alt',
        'N_Reads', 'N_Error', 'N_Low_Quality',
        'N_Inconsistence_Pairs'
    )

    with open(o.query) as fin:
        for line in fin:
            if line.startswith('#'):
                continue

            chr, pos = line.split('\t')[:2]

            reads = bam.get_reads(o, chr, pos)
            unique_pairs, unique_single, nsum, nerror, nlowq, ninconsis = aggregate.aggregate_reads(o, reads)

            d = {'A': 0, 'T': 0, 'C': 0, 'G': 0}

            for x in (unique_single, unique_pairs):
                for *_, base, _ in x:
                    d[base] += 1

            refbase = max(d, key=lambda x: d[x])

            unique_overlapped_pairs_support_alt     = sum(1 for *_, base, paired in unique_pairs if base!=refbase and paired)
            unique_non_overlapped_pairs_support_alt = sum(1 for *_, base, paired in unique_pairs if base!=refbase and not paired)
            unique_single_support_alt               = sum(1 for *_, base, _ in unique_single if base!=refbase)

            println(
                chr, pos,
                unique_overlapped_pairs_support_alt,
                unique_non_overlapped_pairs_support_alt,
                unique_single_support_alt,
                nsum, nerror, nlowq, ninconsis
            )
