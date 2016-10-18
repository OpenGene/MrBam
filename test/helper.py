from random import choice
from os import system as run

def make_bam(dir, s):
    """
    generates test.sam, test.bam, test.bai and test.vcf under $dir by $s

    underline: soft clipped
    dot: base that support ref
    star: base that support alt

    example input:
        r1 + _...........
        r1 - _....*......
        r2 +  ....*......
        r3 +    ......*.....
        r3 -         .*..........
        r4 +    ___...*.....
        r4 -         .*..........
    """

    class Read(object):
        def __init__(self, **args):
            self.__dict__ = args

    def parse_reads():
        for line in s.splitlines():
            for strand in '+', '-':
                if strand in line:
                    name, raw_seq = line.split(strand)
                    seq = raw_seq.lstrip()

                    yield Read(
                        name   = name.strip(),
                        strand = strand,
                        start  = len(raw_seq) - len(seq),
                        end    = len(raw_seq),
                        seq    = list(seq.rstrip())
                    )

                    break

    def generate_random_sequence():
        remainder = { 'A': 'TCG', 'T': 'ACG', 'C': 'ATG', 'G': 'ATC' }
        ref = [ choice('ATCG') for i in range(total_length) ]
        alt = [ choice(remainder[ref[i]]) for i in range(total_length) ]
        return ref, alt

    def calc_cigar():
        for r in reads:
            state, acc, cigar = '\0', 0, ''

            for c in r.seq:
                op = 'S' if c == '_' else 'M'
                if op == state:
                    acc += 1
                else:
                    if acc > 0: # not the initial state
                        cigar += str(acc) + state
                    state, acc = op, 1

            r.cigar = cigar + str(acc) + op

    def calc_pos():
        for r in reads:
            t = ''.join(r.seq).lstrip('_')
            r.pos = r.start + len(r.seq) - len(t)
            r.len = len(t.rstrip('_'))

    def gen_vcf():
        mut = [0] * total_length

        for r in reads:
            for (i,c) in enumerate(r.seq):
                if c == '*':
                    mut[r.start + i] += 1

        with open(dir + "/test.vcf", "w") as fout:
            fout.write("##fileformat=VCFv4.1\n")
            fout.write("##source=MrBam\n")
            fout.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tFORMAT\tNORMAL\tTUMOR\n")

            for i in range(total_length):
                if mut[i] > 0:
                    print(
                        'ref', i, '.', ref[i], alt[i], '.', 'PASS', 'AD', 0, mut[i],
                        file=fout, sep='\t'
                    )

    def replace_mask():
        for i in range(max( r.end for r in reads )):
            for r in reads:
                if r.start <= i < r.end:
                    j = i - r.start
                    r.seq[j] = alt[i] if r.seq[j] == '*' else ref[i]

    def find_mates():
        d = {} # name -> read

        for r in reads:
            if r.name in d: # found pair
                d[r.name].mate = r
                r.mate = d[r.name]
            else:
                d[r.name] = r

    def is_paired(read):
        return 'mate' in read.__dict__

    def calc_tlen():
        for r in reads:
            if is_paired(r):
                if r.strand == '+':
                    r.tlen = r.mate.pos + r.mate.len - r.pos
                else:
                    r.tlen = r.mate.pos - r.pos - r.len

    def calc_flag():
        for r in reads:
            r.flag = 0x03
            r.flag |= 0x40 if r.strand == '+' else 0x80

    def write_sam():
        with open(dir + "/test.sam", 'w') as fout:
            fout.write("@HD\tVN:1.5\tSO:unsorted\n")
            fout.write("@SQ\tSN:ref\tLN:%d\n" % (total_length))

            for r in reads:
                print(
                    r.name, # QNAME
                    r.flag, # FLAG
                    'ref', # RNAME
                    r.pos, # POS
                    60, # MAPQ
                    r.cigar, # CIGAR
                    'ref' if is_paired(r) else '*', # RNEXT
                    r.mate.pos if is_paired(r) else 0, # PNEXT
                    r.tlen if is_paired(r) else 0, # TLEN
                    ''.join(r.seq), # SEQ
                    chr(60 + 33) * len(r.seq), # QUAL
                    file=fout, sep='\t'
                )

    def gen_bam_and_bai():
        run("samtools view -b {0}/test.sam | samtools sort -o {0}/test.bam -".format(dir))
        run("samtools index {0}/test.bam {0}/test.bai".format(dir))

    reads = list(parse_reads())
    total_length = max( r.end for r in reads )
    ref, alt = generate_random_sequence()
    calc_cigar()
    calc_pos()
    gen_vcf()
    replace_mask()
    find_mates()
    calc_tlen()
    calc_flag()
    write_sam()
    gen_bam_and_bai()
