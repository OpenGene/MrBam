from bam import get_reads
from aggregate import aggregate_reads
from enum import Enum

class State(Enum):
    start   = 0
    head    = 1
    comment = 2
    body    = 3

def anno(o):
    "annotate vcf file, adding basic infomations into INFO columns"

    with open(o.query) as fin, open(o.output, 'w') as fout:
        def dispatch(state, line):
            if line.startswith('##'):
                if state in (State.start, State.head):
                    return copy, State.head
                else:
                    raise Exception("unexpected " + line)
            elif line.startswith('#'):
                if state == State.head:
                    return add_head, State.comment
                else:
                    return copy, State.comment
            else:
                return anno, State.body

        def copy(line):
            fout.write(line)

        def add_head(line):
            fout.write("""##FORMAT=<ID=UDP,Number=1,Type=String,Description="Unique_Overlapped_Pairs_Support_Ref, Unique_Non_Overlapped_Pairs_Support_Ref, Unique_Single_Support_Ref, Unique_Overlapped_Pairs_Support_Alt, Unique_Non_Overlapped_Pairs_Support_Alt, Unique_Single_Support_Alt">\n""")
            fout.write(line)

        def anno(line):
            line = line.rstrip().split('\t')

            chr, pos, _, ref, alt = line[:5]
            alt = alt.split('/')

            line[-3] += ':UDP' # format
            line[-2] += ':' # gdna
            line[-1] += ':' # cfdna

            for i, sam in enumerate((o.cfdna, o.gdna)):
                if sam != None:
                    reads = get_reads(o, sam, chr, pos)
                    unique_pairs, unique_single, *_ = aggregate_reads(o, reads)
                    line[-i-1] += ','.join(map(str, count_different_type(unique_pairs, unique_single, alt, ref)))

            print(file=fout, sep='\t', *line)

        state = State.start

        for line in fin:
            action, state = dispatch(state, line)
            action(line)

def count_different_type(unique_pairs, unique_single, alt, ref):
    opr, npr, sr, opa, npa, sa = 0, 0, 0, 0, 0, 0

    for *_, base, paired in unique_pairs:
        if base in alt:
            if paired:
                opa += 1
            else:
                npa += 1
        elif base == ref:
            if paired:
                opr += 1
            else:
                npr += 1

    for *_, base, _ in unique_single:
        if base in alt:
            sa += 1
        elif base == ref:
            sr += 1

    return opr, npr, sr, opa, npa, sa
