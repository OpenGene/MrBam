#!/usr/bin/env python3

from os.path import splitext
from argparse import ArgumentParser, RawDescriptionHelpFormatter
from pysam import AlignmentFile
from MrBam.anno import anno

class SingleMetavarHelpFormatter(RawDescriptionHelpFormatter):
    def _format_action_invocation(self, action):
        if not action.option_strings:
            metavar, = self._metavar_formatter(action, action.dest)(1)
            return metavar
        else:
            parts = []
            if action.nargs == 0:
                parts.extend(action.option_strings)
            else:
                default = action.dest.upper()
                args_string = self._format_args(action, default)
                parts.extend(action.option_strings)
                parts[-1] += ' %s' % args_string
            return ', '.join(parts)

def parse_args():
    description = "example:\n  $ MrBam sample.vcf --cfdna sample_cfdna.bam -o sample_MrBam.vcf --simple"

    parser = ArgumentParser(formatter_class=SingleMetavarHelpFormatter, description=description)

    parser.add_argument('query', help="vcf file contains mutations to query")
    parser.add_argument('-c', '--cfdna', help="bam file contains cfdna reads info. There must be a corresponding .bai file in the same directory")
    parser.add_argument('-g', '--gdna', help="bam file contains gdna reads info. There must be a corresponding .bai file in the same directory")
    parser.add_argument('-o', '--output', help="output vcf file. Will be overwritten if already exists")
    parser.add_argument('-i', '--info', default="/dev/null", help="additional infomations about these position")
    parser.add_argument('-q', '--qual', type=int, default=20, help="drop bases whose qulity is less than this (default: 20)")
    parser.add_argument('-s', '--simple', action='store_true', help="annotate less infomations into vcf output")
    parser.add_argument('-f', '--fast', action='store_true', help="do not infer origin read size by CIGAR, it can be faster and consume less memory.")
    parser.add_argument('--allow-inconsist', action='store_true', help="allow different reads stack at the same position. This increases sensitivity.")
    parser.add_argument('-v', '--verbos', action='store_true', help="output debug info")

    return parser.parse_args()

def init(o):
    if o.cfdna == None and o.gdna == None:
        raise Exception("At least one of --cfdna and --gdna should be specified")

    if o.cfdna != None:
        o.cfdna = AlignmentFile(o.cfdna, "rb")
        if not o.cfdna.has_index():
            raise Exception("Index not found, use `samtools index` to generate")

    if o.gdna != None:
        o.gdna = AlignmentFile(o.gdna, "rb")
        if not o.gdna.has_index():
            raise Exception("Index not found, use `samtools index` to generate")

    if o.output == None:
        basename, extname = splitext(o.query)
        o.output = basename + "_MrBam" + extname

if __name__ == '__main__':
    o = parse_args()
    init(o)
    anno(o)
