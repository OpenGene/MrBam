from MrBam.tools import memo, try_append

def get_reads(o, sam, chr, pos):
    "get all reads covers chr:pos"

    pos = int(pos) - 1 # 0-based leftmost coordinate

    for read in sam.fetch(chr, pos, pos+1):
        aligned_pairs = read.get_aligned_pairs()

        if len(aligned_pairs) == 0:
            if o.verbos:
                print("read not aligned: " + read.query_name)
            continue

        try:
            i, query_pos = next((i, qpos) for (i, (qpos, rpos)) in enumerate(aligned_pairs) if rpos == pos)
            if query_pos == None:
                t = 'D'
            elif i+1 < len(aligned_pairs) and aligned_pairs[i+1][1] == None:
                t = 'I'
            else:
                t = 'M'
        except:
            continue

        yield (
            read.query_name,
            t  if t in ('D', 'I') else read.query_sequence[query_pos],
            -1 if t in ('D', 'I') else read.query_qualities[query_pos],
            read.reference_start - read.query_alignment_start,
            read.infer_query_length(),
            -1 if o.mismatch_limit == -1 else nmismatch(read),
            read.next_reference_start,
            abs(read.template_length),
            read.is_reverse,
            read.is_paired and not read.mate_is_unmapped
        )

def nmismatch(read):
    try:
        md = read.get_tag('MD')
        
        deletion, snp = False, 0
        for i in md:
            if i == '^':
                deletion = True
            elif i in 'ATCG':
                snp += not deletion
            else:
                deletion = False

        return snp + sum(1 for op, length in read.cigartuples if op in (1, 2))

    except KeyError:
        pass

    try:
        nm    = read.get_tag('NM')
        indel = sum(length-1 for op, length in read.cigartuples if op in (1, 2))
        
        return nm - indel

    except KeyError:
        pass

    except TypeError:
        pass

    return -1

@memo
def pad_softclip(sam):
    "return a dict of name -> startpos, length, where softclipped bases were padded"

    namedict, pairdict = {}, {}

    for read in sam.fetch():
        if not read.is_paired or read.mate_is_unmapped:
            continue

        adjusted_start = read.reference_start - read.query_alignment_start
        adjusted_end   = adjusted_start + read.infer_query_length()
        name           = read.query_name

        try_append(namedict, name, (adjusted_start, adjusted_end))

    for k, v in namedict.items():
        if len(v) != 2:
            pairdict[k] = -1, len(v)
        else:
            start  = min(map(lambda x: x[0], v))
            length = max(map(lambda x: x[1], v)) - start
            pairdict[k] = start, length

    return pairdict
