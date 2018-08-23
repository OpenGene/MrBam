from collections import Counter
from MrBam.tools import try_append 
#from tools import try_append

def snv_mut(reads, ref_set, o, pad_softclip = None):
#def snv_mut(reads, ref_set, pad_softclip = None):    
    "gettting all possible mutation combinations, including MNV(Multiple Necluotide Variation)"
    "aggregate reads by startpos, endpos and base"

    snv = []
    name_dict     = {} # name -> reads
    unique_pairs  = {} # start, length, is_overlap -> [base, quality]
    unique_single = {} # start, length, is_reverse -> [base, quality]

    base, mut, qual, start_self, query_len, nmismatch, XA, start_next = ['',''], ['',''], ['',''], ['',''], ['',''], ['',''], ['',''], ['','']
    template_len, reverse, paired, q10, terminal, cigartuples = ['',''], ['',''], ['',''], ['',''], ['',''], ['','']
    reference_start, copy_number, mapping_quality =  ['',''], ['',''], ['','']

    for read in reads:
        name, *info = read
        try_append(name_dict, name, info)

    for name in list(name_dict.keys()):
        if len(name_dict[name]) == 1: # non-overlap or single

            [base[0], mut[0], qual[0], start_self[0], query_len[0], nmismatch[0], XA[0], template_len[0],
             start_next[0], reverse[0], paired[0], q10[0], terminal[0], cigartuples[0], reference_start[0],
             copy_number[0], mapping_quality[0]] = name_dict[name][0]

            if 0 <= qual[0] <= o.qual:
                if o.verbos:
                    print("low quality: " + name)
                del name_dict[name]
                continue

            if o.mismatch_limit != -1:
                if nmismatch[0] > o.mismatch_limit:
                    if o.verbos:
                        print("%s has %d mismatch (limit: %d)" % (name, nmismatch[0], o.mismatch_limit))
                    del name_dict[name]
                    continue

            if o.dropXA and XA[0]:
                if o.verbos:
                    print("multiple alignment: " + name)
                del name_dict[name]
                continue
            
            if o.snp:
                snv.append(''.join(mut[0]))

            if paired[0]:
                if o.fast:
                    start = min(start_self[0], start_next[0])
                else:
                    start, template_len[0] = adjusted_pos[name]
                    if start < 0:
                        if o.verbos:
                            print("%s: more than 2 reads (%d total) share the same name; all droped." % (name, template_len[0]))
                        continue
                if o.snp:
                    try_append(unique_pairs, (start, template_len[0], False), (''.join(mut[0]), qual[0]))
                else:
                    try_append(unique_pairs, (start, template_len[0], False), (base[0], qual[0]))
            else:
                if o.snp:
                    try_append(unique_single, (start_self[0], query_len[0], reverse[0]), (''.join(mut[0]), qual[0]))
                else:
                    try_append(unique_single, (start_self[0], query_len[0], reverse[0]), (base[0], qual[0]))

        elif len(name_dict[name]) == 2:  # pe reads with mutation site in overlap area
            
            [base[0], mut[0], qual[0], start_self[0], query_len[0], nmismatch[0], XA[0], template_len[0],
             start_next[0], reverse[0], paired[0], q10[0], terminal[0], cigartuples[0],reference_start[0],
             copy_number[0], mapping_quality[0]] = name_dict[name][0]

            [base[1], mut[1], qual[1], start_self[1], query_len[1], nmismatch[1], XA[1], template_len[1],
             start_next[1], reverse[1], paired[1], q10[1], terminal[1], cigartuples[1], reference_start[1],
             copy_number[1], mapping_quality[1]]= name_dict[name][1]

            if 0 <= qual[0] <= o.qual and 0 <= qual[1] <= o.qual:
                if o.verbos:
                    print("low quality: " + name)
                del name_dict[name]
                continue

            if base[0] != base[1]:
                if o.verbos:
                    print("pair inconsistent: " + name)
                del name_dict[name]
                continue

            if o.mismatch_limit != -1:
                if nmismatch[0] > o.mismatch_limit or nmismatch[1] > o.mismatch_limit:
                    if o.verbos:
                        print("%s has %d mismatch (limit: %d)" % (name, max(nmismatch[0], nmismatch[1]), o.mismatch_limit))
                    del name_dict[name]                    
                    continue

            if o.dropXA and (XA[0] or XA[1]):
                if o.verbos:
                    print("multiple alignment: " + name)
                del name_dict[name]
                continue                
            
            if o.snp:
                snv.append(''.join(mut[0]))
                snv.append(''.join(mut[1]))

            start = min(start_self[0], start_self[1])
            tlen  = max(start_self[0] + query_len[0], start_self[1] + query_len[1]) - start
            qual[0]  = max(qual[0], qual[1])
            if o.snp:
                try_append(unique_pairs, (start, tlen, True), (''.join(mut[0]), qual[0]))
            else:
                try_append(unique_pairs, (start, tlen, True), (base[0], qual[0]))
        else:
            if o.verbos:
                print("%s: more than 2 reads (%d total) share the same name; all droped." % (name, len(name_dict[name])))

    if o.snp:
        c = Counter(snv)
        snv = {}

        if len(ref_set) == 2:
            # continuous bases
            for mut, num in c.items():
                if mut == ref_set[0] + ref_set[1] or mut == ref_set[0] + 'N' or mut == 'N' + ref_set[1]:
                    continue

                elif mut[0] == 'N':
                    snv['N' + mut[1]] = num

                elif mut[1] == 'N':
                    snv[mut[0] + 'N'] = num

                elif mut[0] != ref_set[0]:

                    if mut[1] == ref_set[1]:
                        snv[mut[0] + 'N'] = num

                    elif mut[1] != ref_set[1]:
                            snv[mut] = num
                    
                elif mut[1] != ref_set[1]:
                        snv['N' + mut[1]] = num
                else:
                    raise Exception ("unclear variation type2 %s" % mut)

        elif len(ref_set) == 1:
            # single variation
            for mut, num in c.items():
                if mut != ref_set[0]:
                    snv[mut] = num
        else:
            raise Exception('len over 2 in ref_set ', ref_set, name_dict)                

        if o.verbos:
            print("initial called snvs: ", snv)
        
        variation = []
        for mut in snv:
            if snv[mut] >= 4:
                # mutation with supporing reads >= 4 will be regarded as novel snv
                variation.append(mut)

        if o.verbos:
            print("initial snvs after filtering ", variation)

        return variation, name_dict, unique_pairs, unique_single

    else:
        return name_dict, unique_pairs, unique_single

