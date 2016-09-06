#!/usr/bin/env python3

from argparse import ArgumentParser
from anno import anno

parser = ArgumentParser()

parser.add_argument('query', help="vcf file contains mutations to query")
parser.add_argument('reads', help="bam file contains origin reads info. There must be a corresponding .bai file in the same directory")
parser.add_argument('--quality', '-q', type=int, default=20, help="drop bases whose qulity is less than this")
parser.add_argument('--ref', '-r', help="reference file. If not provided, consider the mode of reads as the ref")
parser.add_argument('--verbos', '-v', action='store_true', help="output debug info")

o = parser.parse_args()

anno(o)
