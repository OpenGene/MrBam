import re
from collections import Counter

def extra_info(o, name_dict, ref, alt):

    MaxInsertSize = 10000
    nq10, nterminal, nmulti, noverlap_pe, noverlap_se, numi, nCN_1, nCN_2, ntotal, nmapqual, ave_mapqual, \
    totalInsertSize, insertSize, ave_insertSize \
    = [0,0], [0,0], [0,0], [0,0], [0,0], ['NA','NA'], [0,0], [0,0], [0,0], [0,0], [0,0], [0,0], [0,0],[0,0]
    if o.UMI:
        numi = [0,0]

    #number of reads with average quality over 10 and mutation located at 20-80% of read and multiple alignment, and overlap
    UMI_set, Alt_set = [],[]
    # uni_seq -> start, end, may including many sequence's information 
    base, mut, XA, q10, terminal, cigartuples, reference_start,copy_number, mapping_quality \
    = ['',''], ['',''], ['',''], ['',''], ['',''], ['',''], ['',''], ['',''], ['','']

    if o.indel:
        if ref =='-' or len(alt) >= 2:
            alt, ref = 'I', 'ATCG'
        elif alt =='-' or len(ref) >= 2:
            alt, ref = 'D', 'ATCG'

    def umi_mismatch(umi1, umi2):
        nmismatch = 0
        if len(umi1) != len(umi2):
            raise Exception("various UMI seq length: %s %s" %( umi1, umi2))
        for i, x in enumerate(umi1):
            if umi1[i] != umi2[i]:
                nmismatch +=1
        if nmismatch >= 2:
            return False
        else:
            return True    

    for name, reads in name_dict.items():
        if len(reads) == 1:
            try:
                mut[0], qual, r1start, r1len, nmismatch, XA[0], insertSize[0], *_, q10[0], terminal[0], cigartuples[0] \
                , reference_start[0], copy_number[0], mapping_quality[0] = reads[0]
                
            except TypeError:
                print("Variable Issue SE %s %s " % (reads ,reads[0]))
                continue

            i = 0
            if mut[0] in alt:
                if terminal[0] is True:
                    nterminal[i] += 1
                if q10[0] is True:
                    nq10[i] += 1
                if XA[0]:
                    nmulti[i] += 1
                if copy_number[0] >= 2:
                    nCN_2[i] += 1
                else:
                    nCN_1[i] += 1
                noverlap_se[i] += 1
                nmapqual[i] += mapping_quality[0]
                ntotal[i] += 1
                if insertSize[0] <= MaxInsertSize:
                    totalInsertSize[0] += insertSize[0]

                #test
                if o.verbos:
                    print('se insertSize ',name, "\t", insertSize[0])

        elif len(reads) == 2:
            r1, r2 = reads
            try:
                mut[0], qual, r1start, r1len, nmismatch, XA[0], insertSize[0], *_, q10[0], terminal[0], cigartuples[0] \
                , reference_start[0], copy_number[0], mapping_quality[0] = reads[0]
                mut[1], qual, r1start, r1len, nmismatch, XA[1], insertSize[1], *_, q10[1], terminal[1], cigartuples[1] \
                , reference_start[1], copy_number[1], mapping_quality[1] = reads[1]

            except TypeError:
                print("Variable Issue PE %s %s %s" %(reads, r1, r2))
                continue

            j = 0
            for i, base1 in enumerate(mut):
                if base1 in alt:
                    if q10[i] is True:
                        nq10[j] += 1
                    if terminal[i] is True:
                        nterminal[j] += 1
                    if XA[i]:
                        nmulti[j] += 1
                    if copy_number[i] >= 2:
                        nCN_2[j] += 0.5
                    else:
                        nCN_1[j] += 0.5
                    nmapqual[j] += mapping_quality[i]
                    ntotal[j] += 1
                    noverlap_pe[j] += 0.5
                    if insertSize[0] == insertSize[1]:
                        if i == 0 and insertSize[0] <= MaxInsertSize:
                            totalInsertSize[0] += insertSize[0]
                    else:
                        print("pair-end owes various insertSize", name)
                        raise Exception("size: ", insertSize)

                    #test
                    if o.verbos and i == 0:
                        print('pe insertSize ',name, "\t",insertSize[0])

    i = 0
    if ntotal[i] != 0:
        ave_mapqual[i] = int(nmapqual[i] / ntotal[i])                    
        ave_insertSize[i] = int(totalInsertSize[i] / (noverlap_pe[i] + noverlap_se[i]))
        if o.verbos:
            print("totalInsertSize: " ,totalInsertSize[i],"\t",noverlap_se[i] + noverlap_pe[i])
    else:
        ave_mapqual[i] = 0
        

    if o.UMI:
        UMI_set = []
        numi_tmp = 0
        for name, info in name_dict.items():
            umi_judge = re.search(r':UMI_(\w+)_(\w+)', name)
            if umi_judge:
                if info[0][0] in alt:
                    umi_judge = re.search(r':UMI_(\w+)_(\w+)', name)
                    umi_seq_for = umi_judge.group(1) + umi_judge.group(2)
                    umi_seq_rev = umi_judge.group(2) + umi_judge.group(1)
                    UMI_set.append((umi_seq_for,umi_seq_rev))
            else:
                raise Exception("Failing to recognize UMI sequence from sequence's name!: %s " % ( name ))

        for i, umi1 in enumerate(UMI_set):
            if umi1 != '' and  i != len(UMI_set) - 1:
                for j in range(i + 1, len(UMI_set)):
                    if UMI_set[j] != '':
                        if umi_mismatch(umi1[0], UMI_set[j][0]) or umi_mismatch(umi1[0], UMI_set[j][1]) or umi_mismatch(umi1[1], UMI_set[j][0]) or  umi_mismatch(umi1[1], UMI_set[j][1]):
                            UMI_set[j] = ''
       
        for umi1 in UMI_set:
            if umi1 != '':
                numi_tmp += 1
        numi[0] = numi_tmp   

    return nq10[0], nterminal[0], nmulti[0], int(noverlap_pe[0]), noverlap_se[0], numi[0], int(nCN_1[0]), int(nCN_2[0]), ave_mapqual[0], ave_insertSize[0]