# MrBam
[![Build Status](https://travis-ci.org/OpenGene/MrBam.svg?branch=master)](https://travis-ci.org/OpenGene/MrBam)
[![Coverage Status](https://coveralls.io/repos/github/OpenGene/MrBam/badge.svg?branch=master)](https://coveralls.io/github/OpenGene/MrBam?branch=master)

For a given mutation, query its mutated reads from a BAM, merge the reads by positions and give the unique count.

# Prerequisites

- Python 3.4+
- Pysam (`$pip install pysam`)

# Usage

```
$ python -m MrBam.main --help
usage: main.py [-h] [-c CFDNA] [-g GDNA] [-o OUTPUT] [-i INFO] [-q QUAL] [-s]
               [-f] [-v] query

example:
  $ MrBam sample.vcf --cfdna sample_cfdna.bam -o sample_MrBam.vcf --simple

positional arguments:
  query                vcf file contains mutations to query

optional arguments:
  -h, --help            show this help message and exit
  -c, --cfdna CFDNA     bam file contains cfdna reads info. There must be a
                        corresponding .bai file in the same directory
  -g, --gdna GDNA       bam file contains gdna reads info. There must be a
                        corresponding .bai file in the same directory
  -o, --output OUTPUT   output vcf file. Will be overwritten if already exists
  --skip SKIP           skip the first N lines
  -q, --qual QUAL       drop bases whose qulity is less than this (default:
                        25)
  -s, --simple          annotate less infomations into vcf output
  -f, --fast            do not infer origin read size by CIGAR, it can be
                        faster and consume less memory.
  --drop-inconsist      drop different reads stack at the same position. This
                        decreases sensitivity.
  --dropXA              drop reads that has XA tag (multiple alignment)
  -m, --mismatch-limit MISMATCH_LIMIT
                        if set, drop reads that has more mismatches than the
                        limit. requires a 'MD' or a 'NM' tag to be present.
  -v, --verbos          output debug info
```

# Performace

```
#sample  option    bam_size(mb)  vcf_lines  CPU_time(s)  Memory(mb)
Sam3     (default) 194           14978      147          1116
Sam3     --fast    194           14978      129          27
Sam2     (default) 655           33702      500          3162
Sam2     --fast    655           33702      417          28
Sam1     (default) 1620          113066     5952         8377
Sam1     --fast    1620          113066     5785         34
Sam4     (default) 2338          648336     49067        9912
Sam4     --fast    2338          648336     60393        36
```

* CPU_time is user + sys
* Memory may vary accroding to system memory pressure
* Test on Intel(R) Xeon(R) CPU E5-2699 v3 @ 2.30GHz
