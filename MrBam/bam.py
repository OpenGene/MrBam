from MrBam.tools import memo, try_append
#from tools import memo, try_append


def get_reads(o, sam, chr, pos_set, ref_set):
#def get_reads(sam, chr, pos_set, ref_set):
    "get all reads covers chr:pos"
    unique_read = {}
    # Filtering reads with same name except for pair-end reads
 
    for i in range(len(pos_set)):
        pos = int(pos_set[i]) - 1  # 0-based leftmost coordinate
        ref = ref_set[i]

        for read in sam.fetch(chr, pos, pos+1):
            aligned_pairs = read.get_aligned_pairs()
            if len(aligned_pairs) == 0:
                if o.verbos:
                    print("read not aligned: " + read.query_name)
                continue
            
            if read.query_name in unique_read:
                if unique_read[read.query_name][0] == 2: 
                    # Not pair ends reads
                    #print(read.query_name,'2',unique_read[read.query_name])
                    continue
                else:
                    if unique_read[read.query_name][2] != read.reference_start - read.query_alignment_start:
                        unique_read[read.query_name] = [2, '']
                   
                    elif unique_read[read.query_name][1] == i:
                        unique_read[read.query_name] = [2, '']

                    else:
                        #print (read.query_name, '1', "sep", unique_read[read.query_name] )
                        continue              
            else:
                unique_read[read.query_name] = [1, i, read.reference_start - read.query_alignment_start]

            aligned_dict = {}
            tmp = -1
            for base, num in read.cigartuples:
                for x in range(1, num + 1):
                    aligned_dict[x + tmp] = base
                tmp += num
            # correction for snv at the last aligned bases where the rest are softclipped bases

            try:
                j, query_pos = next((j, qpos) for (j, (qpos, rpos)) in enumerate(aligned_pairs) if rpos == pos)
                if query_pos is None:
                    t = 'D'
                elif j+1 < len(aligned_pairs) and aligned_pairs[j+1][1] is None:
                    if aligned_dict[j+1] == 1:
                        # insertion
                        t = 'I'
                    elif aligned_dict[j+1] == 4:
                        # softclipped bases
                        t = 'M'
                    else:
                        print('error in softclipped bases', read.query_name, aligned_pairs, aligned_dict, pos, j)    
                else:
                    #if read.query_sequence[query_pos] != ref:
                    #   t = 'snv'
                    #else:
                        t = 'M'    

            except:
                print('error in location', read.query_name, pos)
                continue

            yield (
                read.query_name,
                t if t in ('D', 'I') else read.query_sequence[query_pos],
                -1 if o.snp is False else continuous_bases(pos_set, read, aligned_pairs, aligned_dict, read.query_sequence, query_pos),
                -1 if t in ('D', 'I') else read.query_qualities[query_pos],
                read.reference_start - read.query_alignment_start,
                read.infer_query_length(),
                #nmismatch(read),
                -1 if o.mismatch_limit == -1 else nmismatch(read),
                read.has_tag("XA"),
                read.next_reference_start,
                abs(read.template_length),
                read.is_reverse,
                read.is_paired and not read.mate_is_unmapped,
                #False if t == 'M' else q10(read),
                #False if t == 'M' else terminal(read, j, t),
                q10(read),
                terminal(read, j, t),
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
    # removing deletions and softclipped bases before mutation site while keeping insertions
    # mutation base wasn't included when calculating 
    if var == "I": 
        x += 1
    loc = x
    dis = 0
    for mut, num in read.cigartuples:
        if mut == 4 or mut == 2:
            for tmp in range(dis, dis + num): # mutation may not be the first base in indel area
                if tmp < x:
                    loc -= 1
                else:
                    loc /= len(read.query_sequence)
                    if 0.1 <= loc < 0.9:
                        return  False
                    else:
                        return  True                       

        elif mut == 0 or mut == 1:
            for tmp in range(dis, dis + num): 
                if tmp == x:
                    loc /= len(read.query_sequence)
                    if 0.1 <= loc < 0.9:
                        return  False
                    else:
                        return  True
        dis += num                         
        
    print("error in mutation location: %s %d " % (read.query_name, x))


def continuous_bases(pos_set, read, aligned_pairs, aligned_dict, read_seq, read_pos):
    if len(pos_set) == 1:
        if read_pos is None:
            return 'D'
        else:
            return read_seq[read_pos]
 
    bases = []
    for pos in pos_set:
        pos = int(pos) - 1
        try:
            x, query_pos = next((x, qpos) for (x, (qpos, rpos)) in enumerate(aligned_pairs) if rpos == pos)
            if query_pos is None:
                bases.append('D')
            elif x+1 < len(aligned_pairs) and aligned_pairs[x+1][1] is None:
                if aligned_dict[x+1] == 1:
                    bases.append('I')
                elif aligned_dict[x+1] == 4:
                    t = read.query_sequence[query_pos]
                    bases.append(t)
            else:
                bases.append(read.query_sequence[aligned_pairs[x][0]])
        except:
            bases.append('N')
            
    return bases                
