def aggregate_reads(o, reads):
    "aggregate reads by startpos, endpos and base"

    name_dict     = {} # name -> reads
    unique_pairs  = {} # start, length, base, is_overlap -> quality
    unique_single = {} # start, length, base, is_reverse -> quality

    nsum,  nerror    = 0, 0 # depth, errors (such as 3 reads share the same name)
    nlowq, ninconsis = 0, 0 # low quality bases, inconsistent pairs (count on reads)

    def try_append(d, k, v):
        if k in d:
            d[k].append(v)
        else:
            d[k] = [v]

    for read in reads:
        name, *info = read
        try_append(name_dict, name, info)
        nsum += 1

    for name, reads in name_dict.items():
        if len(reads) == 1: # non-overlap or single
            base, qual, r1start, r1len, r2start, tlen, isrev, paired = reads[0]

            if qual <= o.quality:
                nlowq += 1
                continue

            tlen = abs(tlen)

            if paired:
                start = min(r1start, r2start)
                try_append(unique_pairs, (start, tlen, base, False), qual)
            else:
                try_append(unique_single, (r1start, tlen, base, isrev), qual)

        elif len(reads) == 2: # overlap
            r1, r2 = reads
            base1, qual1, r1start1, r1len, r2start1, tlen1, isrev1, paired1 = r1
            base2, qual2, r1start2, r1len, r2start2, tlen2, isrev2, paired2 = r2

            if qual1 <= o.quality or qual2 <= o.quality:
                nlowq += 2
                continue

            if base1 != base2:
                ninconsis += 2
                continue

            start = min(r1start1, r2start1)
            tlen  = abs(tlen1)
            qual  = max(qual1, qual2)

            try_append(unique_pairs, (start, tlen, base1, True), qual)

        else: # error
            if o.verbos:
                print("%s: more than 2 reads (%d total) share the same name; all droped." % (name, len(reads)))
            nerror += len(reads)

    return unique_pairs, unique_single, nsum, nerror, nlowq, ninconsis
