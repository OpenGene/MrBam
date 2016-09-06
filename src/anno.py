from os.path import splitext
from bam import get_reads
from aggregate import aggregate_reads

def anno(o):
    "annotate vcf file, adding basic infomations into INFO columns"

    basename, extname = splitext(o.query)

    with open(o.query) as fin, open(basename+"_MrBam"+extname, 'w') as fout, open(basename+"_MrBam.info.txt", 'w') as finfo:
        println = lambda f, *args: f.write('\t'.join(map(str, args)) + '\n')

        println(finfo,
            '#Chr', 'Position',
            'Unique_Overlapped_Pairs_Support_Alt',
            'Unique_Non_Overlapped_Pairs_Support_Alt',
            'Unique_Single_Support_Alt',
            'N_Reads', 'N_Error', 'N_Low_Quality',
            'N_Inconsistence_Pairs'
        )

        if extname != ".txt":
            for line in fin: # copy headings
                if not line.startswith('##'):
                    fout.write("""##FORMAT=<ID=UDP,Number=1,Type=String,Description="Unique_Overlapped_Pairs_Support_Ref, Unique_Non_Overlapped_Pairs_Support_Ref, Unique_Single_Support_Ref, Unique_Overlapped_Pairs_Support_Alt, Unique_Non_Overlapped_Pairs_Support_Alt, Unique_Single_Support_Alt">""")
                    break

                fout.write(line)

        for line in fin:
            line = line.rstrip().split('\t')

            chr, pos, _, ref, alt = line[:5]

            alt = alt.split('/')

            reads = get_reads(o, chr, pos)
            unique_pairs, unique_single, nsum, nerror, nlowq, ninconsis = aggregate_reads(o, reads)

            unique_overlapped_pairs_support_alt     = sum(1 for *_, base, paired in unique_pairs if base in alt and paired)
            unique_overlapped_pairs_support_ref     = sum(1 for *_, base, paired in unique_pairs if base == ref and paired)
            unique_non_overlapped_pairs_support_alt = sum(1 for *_, base, paired in unique_pairs if base in alt and not paired)
            unique_non_overlapped_pairs_support_ref = sum(1 for *_, base, paired in unique_pairs if base == ref and not paired)
            unique_single_support_alt               = sum(1 for *_, base, _ in unique_single if base in alt)
            unique_single_support_ref               = sum(1 for *_, base, _ in unique_single if base == ref)

            line[-3] += ":UDP" # format
            line[-2] += ":"
            line[-1] += ":" + ','.join(map(str, (unique_overlapped_pairs_support_ref,
                                                 unique_non_overlapped_pairs_support_ref,
                                                 unique_single_support_ref,
                                                 unique_overlapped_pairs_support_alt,
                                                 unique_non_overlapped_pairs_support_alt,
                                                 unique_single_support_alt)))

            println(fout, *line)

            println(finfo,
                chr, pos,
                unique_overlapped_pairs_support_alt,
                unique_non_overlapped_pairs_support_alt,
                unique_single_support_alt,
                nsum, nerror, nlowq, ninconsis
            )
