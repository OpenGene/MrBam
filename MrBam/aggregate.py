from MrBam.tools import try_append

def aggregate_reads(o, reads, adjusted_pos=None):
    "aggregate reads by startpos, endpos and base"

    name_dict     = {} # name -> reads
    unique_pairs  = {} # start, length, is_overlap -> [base, quality]
    unique_single = {} # start, length, is_reverse -> [base, quality]

    nsum,  nerror    = 0, 0 # depth, errors (such as 3 reads share the same name)
    nlowq, ninconsis = 0, 0 # low quality bases, inconsistent pairs (count on reads)

    for read in reads:
        name, *info = read
        try_append(name_dict, name, info)
        nsum += 1

    for name, reads in name_dict.items():
        if len(reads) == 1: # non-overlap or single
            base, qual, r1start, r1len, nmismatch, XA, r2start, tlen, isrev, paired = reads[0]

            if 0 <= qual <= o.qual:
                nlowq += 1
                if o.verbos:
                    print("low quality: " + name)
                continue

            if o.mismatch_limit != -1:
                if nmismatch > o.mismatch_limit:
                    nlowq += 1
                    if o.verbos:
                        print("%s has %d mismatch (limit: %d)" % (name, nmismatch, o.mismatch_limit))
                    continue

            if XA:
                nlowq += 1
                if o.verbos:
                    print("multiple alignment: " + name)
                continue

            if paired:
                if o.fast:
                    start = min(r1start, r2start)
                else:
                    start, tlen = adjusted_pos[name]
                    if start < 0:
                        if o.verbos:
                            print("%s: more than 2 reads (%d total) share the same name; all droped." % (name, tlen))
                        continue
                try_append(unique_pairs, (start, tlen, False), (base, qual))
            else:
                try_append(unique_single, (r1start, r1len, isrev), (base, qual))

        elif len(reads) == 2: # overlap
            r1, r2 = reads
            r1base, r1qual, r1start, r1len, r1nmismatch, r1XA, *_ = r1
            r2base, r2qual, r2start, r2len, r2nmismatch, r2XA, *_ = r2

            if 0 <= r1qual <= o.qual and 0 <= r2qual <= o.qual:
                nlowq += 2
                if o.verbos:
                    print("low quality: " + name)
                continue

            if r1base != r2base:
                ninconsis += 2
                if o.verbos:
                    print("pair inconsistent: " + name)
                continue

            if o.mismatch_limit != -1:
                if r1nmismatch > o.mismatch_limit or r2nmismatch > o.mismatch_limit:
                    nlowq += 2
                    if o.verbos:
                        print("%s has %d mismatch (limit: %d)" % (name, nmismatch, o.mismatch_limit))
                    continue

            if r1XA or r2XA:
                nlowq += 2
                if o.verbos:
                    print("multiple alignment: " + name)
                continue

            start = min(r1start, r2start)
            tlen  = max(r1start+r1len, r2start+r2len) - start
            qual  = max(r1qual, r2qual)

            try_append(unique_pairs, (start, tlen, True), (r1base, qual))

        else: # error
            if o.verbos:
                print("%s: more than 2 reads (%d total) share the same name; all droped." % (name, len(reads)))
            nerror += len(reads)

    return unique_pairs, unique_single, nsum, nerror, nlowq, ninconsis
