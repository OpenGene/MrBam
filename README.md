# MrBam
For a given mutation, query its mutated reads from a BAM, merge the reads by positions and give the unique count.

# Prerequisites

- Python 3.4+
- Pysam (`$pip install pysam`)

# Usage

```
$ ./MrBam/src/main.py --help
usage: main.py [-h] [-c CFDNA] [-g GDNA] [-o OUTPUT] [-i INFO] [-q QUAL] [-v]
               query

positional arguments:
  query                vcf file contains mutations to query

optional arguments:
  -h, --help           show this help message and exit
  -c, --cfdna CFDNA    bam file contains cfdna reads info. There must be a
                       corresponding .bai file in the same directory
  -g, --gdna GDNA      bam file contains gdna reads info. There must be a
                       corresponding .bai file in the same directory
  -o, --output OUTPUT  output vcf file. Will be overwrite if already exists
  -i, --info INFO      additional infomations about these position
  -q, --qual QUAL      drop bases whose qulity is less than this (default: 20)
  -v, --verbos         output debug info
```

# Under development
