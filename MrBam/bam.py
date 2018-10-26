from MrBam.tools import memo

def get_infor(o, reads, pos, ref):
    for read in reads:
        aligned_pairs = read.get_aligned_pairs()
        if len(aligned_pairs) == 0:
            if o.verbos:
                print("read not aligned: " + read.query_name)
            continue
      
        aligned_dict = {}
        tmp = -1
        for base, num in read.cigartuples:
            for x in range(1, num + 1):
                aligned_dict[x + tmp] = base
            tmp += num
        # the last matched base with softclipped bases left shouldn't be seen as insertion 

        readfilter = 'false'
        if o.continous:
            for i in range(len(aligned_pairs)):
                if aligned_pairs[i][-1] is not None:
                    if pos[0] < aligned_pairs[i][-1]:
                        #print(pos[0],"filter1",aligned_pairs)
                        readfilter = 'true' 
                    break
            
            for j in range(1,len(aligned_pairs)+1):
                i = -j
                if aligned_pairs[i][-1] is not None:
                    if pos[-1] > aligned_pairs[i][-1]:
                        #print(pos[-1],"filter2",aligned_pairs)
                        readfilter = 'true'
                    break  

        if readfilter == 'false':     
            try:
                j, query_pos = next((j, qpos) for (j, (qpos, rpos)) in enumerate(aligned_pairs) if rpos == pos[0])
                if query_pos is not None:
                    if j + 1 < len(aligned_dict) and aligned_dict[j+1] == 1:
                        t = 'I'
                    elif len(pos) == 1:
                        base = read.query_sequence[query_pos] 
                        qual = read.query_qualities[query_pos] 
                        if base == ref:
                            t = 'M'
                        else:
                            t = 'MIS'
                    else:
                        base = read.query_sequence[query_pos:query_pos+2] 
                        qual = sum(read.query_qualities[query_pos:query_pos+2]) / 2
                        if base == ref:
                            t = 'M'
                        else:
                            t = 'MIS'
                else:
                    t = 'D'

            except:
                raise Exception('error in location', read.query_name, pos, query_pos)

            yield (
                read.query_name,
                t if t in ('D', 'I') else base,
                -1 if t in ('D', 'I') else qual,
                read.reference_start - read.query_alignment_start,
                read.infer_query_length(),
                -1 if o.mismatch_limit == -1 else nmismatch(read),
                read.has_tag("XA"),
                abs(read.template_length),
                read.next_reference_start,
                read.is_reverse,
                read.is_paired and not read.mate_is_unmapped,
                -1 if t == 'M' else q10(read),
                -1 if t == 'M' else terminal(read, j, t),
                read.cigartuples,
                read.reference_start,
                read.get_tag('CN'),
                read.mapping_quality
            )


def nmismatch(read):
    try:
        md = read.get_tag('MD')
        
        deletion, snp = False, 0
        for x in md:
            if x == '^':
                deletion = True
            elif x in 'ATCG':
                snp += not deletion
            else:
                deletion = False

        return snp + sum(1 for op, length in read.cigartuples if op in (1, 2))

    except KeyError:
        pass

    try:
        nm = read.get_tag('NM')
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

    for y, v in namedict.items():
        if len(v) != 2:
            pairdict[y] = -1, len(v)
        else:
            start  = min(map(lambda x: x[0], v))
            length = max(map(lambda x: x[1], v)) - start
            pairdict[y] = start, length

    return pairdict


def q10(read):
    aver = sum(read.query_qualities) / len(read.query_qualities)
    if aver >= 10:
        return True
    else:
        return False


def terminal(read, x, var):
    # removing deletions and softclipped bases before mutation site while keeping insertions left
    # mutation base wasn't included when calculating 
    if var == "I": 
        x += 1
    loc = x + 1
    dis = -1
    sfc = sum(num for mut, num in read.cigartuples if mut == 4 )
    finalLen = len(read.query_sequence) - sfc

    for mut, num in read.cigartuples:
        if mut == 4 or mut == 2:
            for tmp in range(dis + 1, dis + num + 1): # mutation may not be the first base in deletion area
                if tmp < x:
                    loc -= 1
                else:
                    loc -= 1
                    loc /= finalLen
                    if 0.1 < loc <= 0.9:
                        return False
                    else:
                        return True

        elif mut == 0 or mut == 1:
            for tmp in range(dis + 1, dis + num + 1): 
                if tmp == x:
                    loc /= finalLen
                    if 0.1 < loc <= 0.9: 
                        return False
                    else:
                        return True
        dis += num   

    raise Exception("error in mutation location: %s %d " % (read.query_name, x))

def continuous_bases(pos, read, aligned_pairs, aligned_dict):
    bases = []
    for pos in pos:
        pos = int(pos) - 1
        try:
            x, query_pos = next((x, qpos) for (x, (qpos, rpos)) in enumerate(aligned_pairs) if rpos == pos)
            if query_pos is None:
                bases.append('D')

            elif x+1 < len(aligned_pairs) and aligned_pairs[x+1][1] is None:

                if aligned_dict[x+1] == 1:
                    bases.append('I')

                elif aligned_dict[x+1] == 4:
                    # the rest are softclipped not insertion
                   bases.append(read.query_sequence[query_pos])

            else:
                bases.append(read.query_sequence[query_pos])
        except:
            bases.append('N')
            
    return bases                
