def aggregate_reads(o, reads):
    "aggregate reads by startpos, endpos and base"

    pair_dict     = {} # name -> reads
    unique_pairs  = {} # start, length, base -> quality
    unique_single = {} # start, length, base -> quality

    nsum,  nerror    = 0, 0 # depth, errors (such as 3 reads share the same name)
    nlowq, ninconsis = 0, 0 # low quality bases, inconsistent pairs (count on reads)

    # 1. find pairs
    for read in reads:
        name, *info = read

        if name in pair_dict:
            if len(pair_dict[name]) != 1: # 0 or 2
                print(name + ": more than 2 reads share the same name; all droped.")
                nerror += len(pair_dict[name]) + 1
                pair_dict[name] = []
            else:
                pair_dict[name].append(info)
        else:
            pair_dict[name] = [info]

        nsum += 1

    # 2. find unique pair/single s
    for pair in pair_dict.values():
        if len(pair) == 1: # single
            base, qual, r1start, r2start, tlen = pair[0]

            start = min(r1start, r2start)
            tlen = abs(tlen)

            if qual <= o.quality:
                nlowq += 1
                continue

            if (start, tlen, base) in unique_single:
                unique_single[(start, tlen, base)].append(qual)
            else:
                unique_single[(start, tlen, base)] = [qual]
        elif len(pair) == 2: # pair
            r1, r2 = pair
            base1, qual1, r1start1, r2start1, tlen1 = r1
            base2, qual2, r1start2, r2start2, tlen2 = r2

            if qual1 <= o.quality or qual2 <= o.quality:
                nlowq += 2
                continue

            if base1 != base2:
                ninconsis += 2
                continue

            start = min(r1start1, r2start1)
            tlen  = abs(tlen1)
            qual  = max(qual1, qual2)
            base  = base1

            if (start, tlen, base) in unique_pairs:
                unique_pairs[(start, tlen, base)].append(qual)
            else:
                unique_pairs[(start, tlen, base)] = [qual]

    return unique_pairs, unique_single, nsum, nerror, nlowq, ninconsis
