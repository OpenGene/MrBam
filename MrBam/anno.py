from MrBam.bam import get_reads, pad_softclip
from MrBam.count import count_different_type
from MrBam.extract import extra_info
from MrBam.snvclassify import snv_mut
import re
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
    #else:
        #print('No repeat file is found!')    

    def repeat_area(rep_ref, chr, pos, o):
        if o.repeat is None:
            return None
        pos = int(pos)
        if chr in rep_ref:
            if o.indel: # indel 
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
                if o.indel:
                    return anno_indel, State.body
                else:
                    return anno_snv, State.body

        def copy(line):
            fout.write(line)

        def add_head(line):
            if o.simple:
                fout.write("""##FORMAT=<ID=UDP,Number=1,Type=String,Description="Unique DNA Supports Alt, Non-Overlaps Are Treated As Single. Comma Separated: Multiple_Overlap, Multiple_Single, One_Overlap, One_Single">\n""")
            else:
                fout.write("""##FORMAT=<ID=UDP,Number=1,Type=String,Description="Unique DNA. Comma Separated: Multiple_Overlap_Ref, Multiple_Nonoverlap_Ref, Multiple_Single_Ref, One_Overlap_Ref, One_Nonoverlap_Ref, One_Single_Ref, Multiple_Overlap_Alt, Multiple_Nonoverlap_Alt, Multiple_Single_Alt, One_Overlap_Alt, One_Nonoverlap_Alt, One_Single_Alt">\n""")
            fout.write(line)
       
        def anno_snv(line_set):
            
            if o.verbos:
                print('line_set ', line_set)
            pos_set, ref_set, ori_alt_set, both_set, chr = [], [], [], [], ''
            for j, line in enumerate(line_set):
                chr = line[0]
                pos_set.append(line[1])
                ref_set.append(line[3])
                ori_alt_set.append(line[4])

            for i, sam in enumerate((o.gdna, o.cfdna)):
                if sam is not None:
                    reads = get_reads(o, sam, chr, pos_set, ref_set)
                    #variation, name_dict, unique_pairs, unique_single= snv_mut(reads, ref_set)
                    both_set.append(snv_mut(reads, ref_set, o, None if o.fast else pad_softclip(sam)))
                   
            all_mut = list(set(both_set[0][0]).union(set(both_set[1][0])))
            for n, snv in enumerate(all_mut):
                # Deletion is not included when calling snp
                if snv == 'D':
                    del all_mut[n]

            repeat = repeat_area(rep_ref, chr, pos_set[0], o)
            if len(all_mut) == 0:
                if len(line_set) == 1:
                    all_mut.append(ori_alt_set[0])
                if len(line_set) == 2:
                    all_mut.append(''.join((ori_alt_set[0], 'N')))
                    all_mut.append(''.join(('N', ori_alt_set[1])))
            
            if o.verbos:
                print('Final merged snv:', all_mut)

            for l, alt in enumerate(all_mut):
                if 'N' not in alt:
                    ref = [''.join(ref_set)]
                    alt_set = [alt]
                    line = line_set[0]

                elif alt.startswith('N'):
                    ref = ('N' + ref_set[1], ''.join(ref_set))
                    alt_set = [alt, ref_set[0] + alt[1]]
                    line = line_set[1]

                elif alt.endswith('N'):
                    ref = (ref_set[0] + 'N', ''.join(ref_set))
                    alt_set =[alt, alt[0] + ref_set[1]]
                    line = line_set[0]
            
                elif re.search(r'^[ATGC]$', alt):
                    alt_set = [alt]
                    ref = ref_set
                    if len(line_set) == 1:
                        line = line_set[0]
                    else:
                        line = line_set[l]
                else:
                    raise Exception("mutation type error: ", alt)
                    
                line[-3] += ':UDP'  # format
                line[-2] += ':'  # gdna
                line[-1] += ':'  # cfdna

                tmp = ['','','','']
                for i in range(2):
                    # gdna cfdna
                    mor, mnr, msr, oor, onr, osr, moa, mna, msa, ooa, ona, osa, _ = count_different_type(o, both_set[i][2], both_set[i][3], alt_set, ref)
                    if o.simple:
                        tmp[i] = line[-2 + i] + ','.join(map(str, (moa, mna + msa, ooa, ona + osa)))
                    else:
                        tmp[i] = line[-2 + i] + ','.join(map(str, (mor, mnr, msr, oor, onr, osr, moa, mna, msa, ooa, ona, osa)))

                    (nq10_a, nterminal_a, nmulti_a, noverlap_pe_a, noverlap_se_a, numi_a, nCN_1_a, nCN_2_a, ave_mapqual_a,
                    nq10_r, nterminal_r, nmulti_r, noverlap_pe_r, noverlap_se_r, numi_r, nCN_1_r, nCN_2_r, ave_mapqual_r) \
                    =  extra_info(alt_set, ref, both_set[i][1], o)
                    extra_1 = ','.join(map(str, (nq10_r, nterminal_r, nmulti_r, noverlap_pe_r, noverlap_se_r, numi_r, nCN_1_r, nCN_2_r, ave_mapqual_r )))
                    extra_2 = ','.join(map(str, (nq10_a, nterminal_a, nmulti_a, noverlap_pe_a, noverlap_se_a, numi_a, nCN_1_a, nCN_2_a, ave_mapqual_a)))
                    tmp[2 + i] = ':'.join((str(repeat), extra_1, extra_2))
                    
                line[3] = re.sub(r'N', '', ref[0])
                if alt != '':
                    line[4] = re.sub(r'N', '', alt)

                print(file=fout, sep='\t', *line[:-2], *tmp)  

        def anno_indel(line):
            line = line.rstrip().split('\t')

            chr, pos, _, ref, alt = line[:5]
            pos, ref, = [pos], [ref]
            tmp = ['', '']

            line[-3] += ':UDP' # format
            line[-2] += ':' # gdna
            line[-1] += ':' # cfdna

            repeat = repeat_area(rep_ref, chr, pos[0], o)
            for i, sam in enumerate((o.cfdna, o.gdna)):
                if sam is not None:
                    reads = get_reads(o, sam, chr, pos, ref)
                    name_dict, unique_pairs, unique_single = snv_mut(reads, ref, o,  None if o.fast else pad_softclip(sam))
                    mor, mnr, msr, oor, onr, osr, moa, mna, msa, ooa, ona, osa, _ = count_different_type(o, unique_pairs, unique_single, alt, ref)
                    if o.simple:
                        line[-i-1] += ','.join(map(str, (moa, mna + msa, ooa, ona + osa)))
                    else:
                        line[-i-1] += ','.join(map(str, (mor, mnr, msr, oor, onr, osr, moa, mna, msa, ooa, ona, osa)))

                    (nq10_a, nterminal_a, nmulti_a, noverlap_pe_a, noverlap_se_a, numi_a, nCN_1_a, nCN_2_a, ave_mapqual_a,
                    nq10_r, nterminal_r, nmulti_r, noverlap_pe_r, noverlap_se_r, numi_r, nCN_1_r, nCN_2_r, ave_mapqual_r) \
                    =  extra_info(alt, ref, name_dict, o)

                    extra_1 = ','.join(map(str, (nq10_r, nterminal_r, nmulti_r, noverlap_pe_r, noverlap_se_r, numi_r, nCN_1_r, nCN_2_r, ave_mapqual_r )))
                    extra_2 = ','.join(map(str, (nq10_a, nterminal_a, nmulti_a, noverlap_pe_a, noverlap_se_a, numi_a, nCN_1_a, nCN_2_a, ave_mapqual_a)))
                    tmp[i] = ':'.join((str(repeat), extra_1, extra_2))

            line.append(tmp[1])  # gdna
            line.append(tmp[0])  #cfdna
            print(file=fout, sep='\t', *line)
              
        state = State.before
        
        if o.snp:   
            ###### merge continuous snv,but not indel
            line_set = []
            line_last = ''
            for line in fin:    
                if re.match(r'chr(\d{1,2}|[XY])\t(\d+)', line):
                    # Mutation site part 
                    now_split = line.rstrip().split('\t')

                    if len(line_set) == 0: 
                        line_set.append(now_split)

                    elif len(line_set) == 1:
                        if now_split[0] == line_set[0][0]:
                             # same chromosome id
                            if int(now_split[1]) < int(line_set[0][1]):
                                raise Exception ("%s was not sorted according to chromosome id and mutation position, try again after sorting vcf file\n" % (query))
                            elif int(now_split[1]) - int(line_set[0][1]) == 1:
                                line_set.append(now_split)
                                anno_snv(line_set)
                                line_set = []

                            else: 
                                anno_snv(line_set)
                                line_set[0] = now_split 
                        else:
                            anno_snv(line_set)
                            line_set[0] = now_split
                    else:
                        raise Exception('error in line_set', line_set)

                else: # header part 
                    action, state = dispatch(state, line)
                    action(line)
            
            if len(line_set) != 0:
                anno_snv(line_set)
            
        elif o.indel:
            for line in fin:
                action, state = dispatch(state, line)
                action(line)

