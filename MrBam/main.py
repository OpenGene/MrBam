#!/usr/bin/env python3

from os.path import splitext
from argparse import ArgumentParser, RawDescriptionHelpFormatter
from pysam import AlignmentFile
from MrBam.anno import anno
from datetime import datetime
from MrBam.continous import continous
from MrBam.explain import explain
import logging

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
    parser.add_argument('--skip', type=int, default=0, help="skip the first N lines")
    parser.add_argument('-q', '--qual', type=int, default=25, help="drop bases whose qulity is less than this (default: 25)")
    parser.add_argument('-s', '--simple', action='store_true', help="annotate less infomations into vcf output")
    parser.add_argument('-f', '--fast', action='store_true', help="do not infer origin read size by CIGAR, it can be faster and consume less memory.")
    parser.add_argument('--drop-inconsist', action='store_true', help="drop different reads stack at the same position. This decreases sensitivity.")
    parser.add_argument('--dropXA', action='store_true', help="drop reads that has XA tag (multiple alignment)")
    parser.add_argument('-m', '--mismatch-limit', type=int, default=-1, help="if set, drop reads that has more mismatches than the limit. requires a 'MD' or a 'NM' tag to be present.")
    parser.add_argument('-v', '--verbos', action='store_true', help="output debug info")
    parser.add_argument('-r', '--repeat', help="repeat region in huam genome")
    parser.add_argument('-u', '--UMI', action='store_true', help="count umi sequences when sample sequenced by duplex UMI")
    parser.add_argument('--indel', action='store_true',help="only indel exists in vcf file")
    parser.add_argument('--snp', action='store_true',help="only snp exists in vcf file")
    parser.add_argument('--alt', action='store_true',help="only count reads'info with snv or indel")
    parser.add_argument('--continous', action='store_true',help="count for continuous mutation site")
    parser.add_argument('--explain', action='store_true', help="detail explanation for result")

    return parser.parse_args()

def init(o):
    if o.explain:
        explain()

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

    if o.indel is False and o.snp is False:
        raise Exception("Updated MrBam requires explicit type of mutation, one of --snp or --indel must be provided")

if __name__ == '__main__':
    t1 = datetime.now()
    o = parse_args()
    init(o)
    if o.continous:
        continous(o)
    else:
        anno(o)
    t2 = datetime.now()
    t_used = (t2 - t1).seconds
    #logging.warning("analysis of %s was finished, %d seconds used !" % (o.query, t_used))
    
