from MrBam.bam import get_infor
from MrBam.count import count_different_type
from MrBam.extract import extra_info
from MrBam.snvclassify import snv_mut
from MrBam.tools import try_append
from copy import deepcopy

def continous(o):
    info = {}
    fout = open(o.output,'w')
    with open(o.query) as fin:
        print(fin.readline(),end ='')
        # First mutation site
        line = fin.readline().rstrip()
        sep = line.split('\t')
        Lastpos = int(sep[1]) -1
        LastChr = sep[0]
        info[Lastpos] = line  
        poss = [Lastpos]

        for line in fin:
            sep = line.rstrip().split('\t')
            chr = sep[0]
            pos = int(sep[1]) - 1
           
            if chr != LastChr:
                if len(poss) >= 2:
                    output(poss,info,o,fout)
                elif len(poss) == 1:
                    Lastpos = pos
                    print(info[poss[0]])
                else:
                    raise Exception(len(poss))
                Lastpos = pos
                LastChr = chr
                poss = [pos]
                info = {}
                info[pos] = line.rstrip() 
                continue       

            if pos - Lastpos == 1:
                poss.append(pos)

            else:
                if len(poss) >= 2:
                    output(poss,info,o,fout)
                elif len(poss) == 1:
                    Lastpos = pos
                    print(info[poss[0]])
                else:
                    raise Exception(len(poss))
                info = {}
                poss = [pos]
            
            Lastpos = pos
            LastChr = chr
            info[pos] = line.rstrip()

        if len(poss) > 1:
            output(poss,info,o,fout)
        elif len(poss) == 1:
            print(info[poss[0]])


def output(poss, info, o, fout):
    if o.verbos:
        print("continus", poss)
    if len(poss) > 2:
        for i in poss:
            print(info[i])
        return 
    sep1 = info[poss[0]].split('\t')
    sep2 = info[poss[1]].split('\t')
    ref = sep1[3] + sep2[3]
    alt = sep1[4] + sep2[4]

    fetch = o.cfdna.fetch(sep1[0],poss[0],poss[1]+1)
    reads = []
    for read in fetch:
        reads.append(read)
    
    new_reads = get_infor(o, reads, poss, ref)
    mut_set, name_dict, unique_pairs, unique_single= snv_mut(new_reads, o, ref, alt, None)
    if len(mut_set) == 0:
        if o.verbos:
            print('continus fail ', poss)
        for i in poss:
            print(info[i])
        return 
        
    sep1[-3] += ':UDP'
    sep1[3] = ref
    sep1[4] = alt
    sep1[2] = int(sep1[1]) + 1

    if o.verbos:
        print("test\t", mut_set, sep1[0], poss[0], ref, sep1[3], sep1[4], alt, len(name_dict),len(unique_pairs),len(unique_single))

    mor, mnr, msr, oor, onr, osr, moa, mna, msa, ooa, ona, osa, _ = count_different_type(o, unique_pairs, unique_single, alt, ref)
    (nq10_a, nterminal_a, nmulti_a, noverlap_pe_a, noverlap_se_a, numi_a, nCN_1_a, nCN_2_a, nNonterminalMolecule_a, \
        ave_mapqual_a, ave_insertSize_a) = extra_info(o, name_dict, ref, alt)
    extra_1 = ','.join(map(str, (nq10_a, nterminal_a, nmulti_a, noverlap_pe_a, noverlap_se_a, numi_a, nCN_1_a, nCN_2_a, nNonterminalMolecule_a, ave_mapqual_a, ave_insertSize_a)))
    sep1[-1] += ':' + ','.join(map(str, (mor, mnr, msr, oor, onr, osr, moa, mna, msa, ooa, ona, osa)))

    fetch = o.gdna.fetch(sep1[0],poss[0],poss[1]+1)
    reads = []
    for read in fetch:
        reads.append(read)
    
    new_reads = get_infor(o, reads, poss, ref)
    mut_set, name_dict, unique_pairs, unique_single= snv_mut(new_reads, o, ref, alt, None)
    mor, mnr, msr, oor, onr, osr, moa, mna, msa, ooa, ona, osa, _ = count_different_type(o, unique_pairs, unique_single, alt, ref)
    (nq10_a, nterminal_a, nmulti_a, noverlap_pe_a, noverlap_se_a, numi_a, nCN_1_a, nCN_2_a, nNonterminalMolecule_a, \
        ave_mapqual_a, ave_insertSize_a) = extra_info(o, name_dict, ref, alt)
    extra_2 = ','.join(map(str, (nq10_a, nterminal_a, nmulti_a, noverlap_pe_a, noverlap_se_a, numi_a, nCN_1_a, nCN_2_a, \
     nNonterminalMolecule_a, ave_mapqual_a, ave_insertSize_a)))
    sep1[-2] += ':' + ','.join(map(str, (mor, mnr, msr, oor, onr, osr, moa, mna, msa, ooa, ona, osa)))

    sep1.append(extra_2)
    sep1.append(extra_1)
    print(file = fout, sep = "\t", *sep1)