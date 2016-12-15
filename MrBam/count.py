from collections import Counter

def count_different_type(o, pairs, single, alt, ref):
    mor, mnr, msr, oor, onr, osr, moa, mna, msa, ooa, ona, osa, inconsis = 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0

    if ref=='-' or any(len(x) > 1 for x in alt):
        alt, ref = 'I', 'ATCG'
    elif alt=='-' or len(ref) >= 2:
        alt, ref = 'D', 'ATCG'

    for (*_, overlaped), reads in pairs.items():
        c = Counter( base for base, qual in reads )

        if not o.drop_inconsist:
            c = c.items()
        else:
            num_all = sum(c.values())
            c = c.most_common(1)
            base, num_major = c[0]

            if num_major / num_all < 0.9:
                if o.verbos:
                    print("reads have different base at the same pos (%d %s, %d others): " % (num_major, base, num_all - num_major))
                inconsis += num_all
                continue

        for base, num in c:
            if num == 1:
                if base in alt:
                    if overlaped:
                        ooa += 1
                    else:
                        ona += 1
                elif base in ref:
                    if overlaped:
                        oor += 1
                    else:
                        onr += 1
            else:
                if base in alt:
                    if overlaped:
                        moa += 1
                    else:
                        mna += 1
                elif base in ref:
                    if overlaped:
                        mor += 1
                    else:
                        mnr += 1

    for reads in single.values():
        c = Counter( base for base, qual in reads )

        if not o.drop_inconsist:
            c = c.items()
        else:
            num_all = sum(c.values())
            c = c.most_common(1)
            base, num_major = c[0]

            if num_major / num_all < 0.9:
                if o.verbos:
                    print("reads have different base at the same pos (%d %s, %d others): " % (num_major, base, num_all - num_major))
                inconsis += num_all
                continue

        for base, num in c:
            if num == 1:
                if base in alt:
                    osa += 1
                elif base in ref:
                    osr += 1
            else:
                if base in alt:
                    msa += 1
                elif base in ref:
                    msr += 1

    return mor, mnr, msr, oor, onr, osr, moa, mna, msa, ooa, ona, osa, inconsis
