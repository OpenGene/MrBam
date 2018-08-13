from MrBam.bam import get_reads, pad_softclip
from MrBam.aggregate import aggregate_reads
from MrBam.count import count_different_type
from enum import Enum

class State(Enum):
    before = 0 # info before ##FORMAT
    under  = 1 # ##FORMAT
    after  = 2 # info after ##FORMAT
    body   = 3

def anno(o):
    "annotate vcf file, adding basic infomations into INFO columns"
    rep_ref = {}
    if o.repeat is not None: 
        ## adding human genome's repeat region reference, mutation around repeat region is usually suspectable
        with open(o.repeat) as rep:
            rep_line = rep.rstrip().split('\t')
            rep_ref[rep_line[0]] = (rep_line[1], rep_line[2])
    else:
        print('No repeat file is found!')    

    def repeat_area(rep_ref, chr, pos, ref, alt):
        pos = int(pos)
        if chr in rep_ref:
            if (ref == '-' or len(alt) > 1) or (alt == '-' or len(ref) >= 2): # indel mutation
                if pos in range(int(rep_ref[chr][0]) -1, int(rep_ref[chr][1]) + 1):
                    return True
                else:
                    return False                             
            else: # snv
                if pos in range(int(rep_ref[chr][0]), int(rep_ref[chr][1]) + 2):
                    return True
                else:
                    return False
        return False

    with open(o.query) as fin, open(o.output, 'w') as fout:
        def dispatch(state, line):
            if o.skip > 0:
                o.skip -= 1
                return copy, State.before
            if line.startswith('##FORMAT'):
                if state in (State.before, State.under):
                    return copy, State.under
                else:
                    raise Exception("unexpected " + line)
            elif line.startswith('#'):
                if state in (State.before, State.after):
                    return copy, state
                elif state == State.under:
                    return add_head, State.after
                else:
                    raise Exception("unexpected " + line)
            else:
                return anno, State.body

        def copy(line):
            fout.write(line)

        def add_head(line):
            if o.simple:
                fout.write("""##FORMAT=<ID=UDP,Number=1,Type=String,Description="Unique DNA Supports Alt, Non-Overlaps Are Treated As Single. Comma Separated: Multiple_Overlap, Multiple_Single, One_Overlap, One_Single">\n""")
            else:
                fout.write("""##FORMAT=<ID=UDP,Number=1,Type=String,Description="Unique DNA. Comma Separated: Multiple_Overlap_Ref, Multiple_Nonoverlap_Ref, Multiple_Single_Ref, One_Overlap_Ref, One_Nonoverlap_Ref, One_Single_Ref, Multiple_Overlap_Alt, Multiple_Nonoverlap_Alt, Multiple_Single_Alt, One_Overlap_Alt, One_Nonoverlap_Alt, One_Single_Alt">\n""")
            fout.write(line)

        def anno(line):
            line = line.rstrip().split('\t')

            chr, pos, _, ref, alt = line[:5]

            line[-3] += ':UDP' # format
            line[-2] += ':' # gdna
            line[-1] += ':' # cfdna 

            for i, sam in enumerate((o.cfdna, o.gdna)):
                if sam != None:
                    reads = get_reads(o, sam, chr, pos, ref)
                    unique_pairs, unique_single, *_, name_dict = aggregate_reads(o, reads, None if o.fast else pad_softclip(sam))
                    mor, mnr, msr, oor, onr, osr, moa, mna, msa, ooa, ona, osa, _ = count_different_type(o, unique_pairs, unique_single, alt, ref)
                    if o.simple:
                        line[-i-1] += ','.join(map(str, (moa, mna + msa, ooa, ona + osa)))
                    else:
                        line[-i-1] += ','.join(map(str, (mor, mnr, msr, oor, onr, osr, moa, mna, msa, ooa, ona, osa)))

                    nq10, nmiddle, nmulti, noverlap_pe, noverlap_se, numi, nCN_1, nCN_2 = extra_info(alt, ref, name_dict, o.umi)

                    repeat = repeat_area(rep_ref, chr, pos, ref, alt)
                    line.append(','.join(map(str, (repeat, nq10, nmiddle, nmulti, noverlap_pe, noverlap_se, numi, nCN_1, nCN_2))))
                 
            print(file=fout, sep='\t', *line)

        state = State.before

        for line in fin:
            action, state = dispatch(state, line)
            action(line)
