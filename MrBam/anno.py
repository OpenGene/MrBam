from MrBam.bam import get_infor
from MrBam.count import count_different_type
from MrBam.extract import extra_info
from MrBam.snvclassify import snv_mut
from MrBam.tools import try_append
from copy import deepcopy

def anno(o):
    "annotate vcf file, adding basic infomations into INFO columns"

    rep_ref = {}
    if o.repeat is not None: 
        ## adding human genome's repeat region reference, mutation around repeat region is usually suspectable
        with open(o.repeat) as rep:
            rep_line = rep.rstrip().split('\t')
            rep_ref[rep_line[0]] = (rep_line[1], rep_line[2])

    mutation = {}
    variation_info = {}

    LastChr = ''
    LastPos = 0

    with open(o.query) as fin:
        # vcf file
        for line in fin:
            line = line.rstrip().split('\t')
            chr = line[0]
            try:
                pos =int(line[1]) - 1
            except:
                raise Exception("vcf header should be expluded", o.query)

            if chr != LastChr:
                LastPos = pos
                LastChr = chr

            else:
                if pos < LastPos:
                    raise Exception("vcf ", o.query, " is unsorted!\n", 'chr: ', chr, "pos: ", LastPos, pos)
                LastPos = pos

            try_append(mutation, chr, (pos,line[4]))
            tmp = ''.join((chr,str(pos),line[4]))
            # various mutation type may happen at one site, adding alteration to recogenize
            variation_info[tmp] = line

    with open(o.output, 'w') as fout:
        for k, sam in enumerate((o.gdna, o.cfdna)):

            variation = deepcopy(mutation) # copy from vcf file
            reads = {} # pos-> read from bam file

            for chr in mutation.keys():
                for pos in mutation[chr]:
                    ID=''.join((chr,str(pos[0]),pos[1]))
                    reads[ID] = []

            LastChr = ''
            LastPos = 0
            for read in sam:

                try:
                    cigarAvalable = read.reference_end - read.reference_start
                except:
                    continue

                chr = read.reference_name
                
                if chr != LastChr:
                    LastPos = read.reference_start
                    LastChr = chr
                else:
                    if read.reference_start < LastPos:
                        raise Exception("sam ", sam, " is unsorted!\n")
                    LastPos = read.reference_start

                if chr in variation:
                        pos = variation[chr][0][0]
                        alt = variation[chr][0][1]
                        ID=''.join((chr,str(pos),alt))
                        
                        if read.reference_start <= pos < read.reference_end:
                        # read aligned with first variation site
                            reads[ID].append(read)

                            try:
                                for j in range(1, len(variation[chr])):
                                    if variation[chr][j][0] < read.reference_end:

                                        pos = variation[chr][j][0]
                                        alt = variation[chr][j][1]
                                        ID2 = ''.join((chr,str(pos),alt))
                                        reads[ID2].append(read)
                                    
                                    else:
                                        break
                            except:
                                raise Exception("error in try1\n")

                        elif read.reference_end <= pos:
                            # -- |
                            continue

                        else:
                            # | --
                            variation_info= output(chr, pos, alt, variation_info, o, fout, reads[ID], k, rep_ref)
                            if len(variation[chr]) == 1:
                                del variation[chr]
                                continue
                            else:
                                del variation[chr][0]
                            reads[ID] = []

                            # Reads should be shorter than 200 bases
                            try:
                                for i in range(0,len(variation[chr])) :
                                    # removing variations not aligned with read
                                    pos = variation[chr][0][0]
                                    alt = variation[chr][0][1]
                                    ID = ''.join((chr,str(pos),alt))

                                    if read.reference_start <= pos < read.reference_end:
                                        reads[ID].append(read) 
                                        # the second variation may aligned with vcf file

                                        for j in range(1, len(variation[chr])):
                                            if variation[chr][j][0] < read.reference_end:
                                                pos = variation[chr][j][0]
                                                alt = variation[chr][j][1]
                                                ID2 = ''.join((chr,str(pos),alt))
                                                reads[ID2].append(read)

                                            else:
                                                break
                                        break

                                    elif read.reference_end <= pos:
                                        # -- |
                                        break

                                    else: # | ---
                                        variation_info= output(chr, pos, alt, variation_info, o, fout, reads[ID], k, rep_ref)
                                        reads[ID] = []
                                        if len(variation[chr]) == 1:
                                            del variation[chr]
                                            break
                                        else:
                                            del variation[chr][0]
                            except:
                                raise Exception("variation ",variation[chr])

            # The last read in one chr aligned with mutation site, leaving mutation info undealed
            for chr in variation:
                for pos in variation[chr]:
                    ID=''.join((chr,str(pos[0]),pos[1]))
                    if reads[ID] != []:
                        variation_info= output(chr, pos, variation_info, o, fout, reads[ID], k, rep_ref)
                    else:
                        raise Exception("mutation more than 1",chr, pos, variation)


def output(chr, pos, alt, variation_info, o, fout, reads_pos, k, rep_ref):
    tmp = ''.join((chr,str(pos), alt))
    if k == 0:
        variation_info[tmp][-3] += ':UDP'
    line = deepcopy(variation_info[tmp])

    ref = variation_info[tmp][3]

    repeat = repeat_area(rep_ref, chr, pos + 1, o)
    pos_set = [pos]
    # to corporate with continous.py
    new_reads = get_infor(o, reads_pos, pos_set, ref)
    mut_set, name_dict, unique_pairs, unique_single= snv_mut(new_reads, o, ref, alt, None if o.fast else pad_softclip(sam))

    for each_mut in mut_set:
        tmp =''.join((chr,str(pos), each_mut))
        if each_mut != alt:
            # MNV(multiple nucleotide variation)
            if k == 0:
                # novel mutation in gdna
                variation_info[tmp] = deepcopy(line)
                variation_info[tmp][4] = each_mut

            else:
                if tmp not in variation_info:
                # cfdna owing novel mutation while gnda not
                # gdna owing novel muatation while cfdna not won't be output
                    variation_info[tmp] = deepcopy(line)
                    #This is  copy from gdna-produced result
                    variation_info[tmp][4] = each_mut

                    try:
                        gdna_alt = variation_info[tmp][-3].split(':')[-1].split(',')
                        for i in range(-6,0):
                            gdna_alt[i] = '0'
                        variation_info[tmp][-3] = ':'.join((variation_info[tmp][-3].split(':')[:-1])) + ":" + ','.join(gdna_alt)
                    except:
                        raise Exception(gdna_alt, i, variation_info[tmp])

                    #extra_2 = repeat + ',' +  ','.join(map(str, [0] * 9))
                    extra_2 =','.join(map(str, [0] * 9))
                    variation_info[tmp][-1] = extra_2

        mor, mnr, msr, oor, onr, osr, moa, mna, msa, ooa, ona, osa, _ = count_different_type(o, unique_pairs, unique_single, each_mut, ref)
        (nq10_a, nterminal_a, nmulti_a, noverlap_pe_a, noverlap_se_a, numi_a, nCN_1_a, nCN_2_a, nNonterminalMolecule_a, ave_mapqual_a, ave_insertSize_a) = extra_info(o, name_dict, ref, each_mut)
        variation_info[tmp][-2] += ':' + ','.join(map(str, (mor, mnr, msr, oor, onr, osr, moa, mna, msa, ooa, ona, osa)))
        extra_2 = ','.join(map(str, (nq10_a, nterminal_a, nmulti_a, noverlap_pe_a, noverlap_se_a, numi_a, nCN_1_a, nCN_2_a, nNonterminalMolecule_a, ave_mapqual_a, ave_insertSize_a)))
        variation_info[tmp].append(extra_2)

        if k == 1:
            print(file=fout, sep='\t', *variation_info[tmp])
            fout.flush()
            del variation_info[tmp]
        
    return variation_info

def repeat_area(rep_ref, chr, pos, o):
    if o.repeat is None:
        return 'NA'

    if chr in rep_ref:
        if int(rep_ref[chr][0]) <= pos < int(rep_ref[chr][1]) + 1:
            return True
        else:
            return False                             
    return False