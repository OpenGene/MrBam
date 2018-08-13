import re
from collections import Counter

def extra_info(alt, ref, name_dict, umi):

    nq10, nmiddle, nmulti, noverlap_pe, noverlap_se, numi, nCN_1, nCN_2 = 0, 0, 0, 0, 0, 'None', 0, 0
    #number of reads with average quality over 10 and mutation located at 20-80% of read and multiple alignment, and overlap
    UMI_set, Alt_set = [],[]
    nerror = 0
    # uni_seq -> start, end, may including many sequence's information 
    #base, XA, q10, middle, cigartuples, reference_start = [None,None], [None,None], [None,None], [None,None], [None,None], [None,None]
    base, XA, q10, middle, cigartuples, reference_start, copy_number = ['',''], ['',''], ['',''], ['',''], ['',''], ['',''], ['','']


    if ref == '-' or len(alt) >= 2:
        alt, ref = 'I', 'ATCG'
    elif alt == '-' or len(ref) >= 2:
        alt, ref = 'D', 'ATCG'

    def umi_mismatch(umi1, umi2):
        nmismatch = 0
        if len(umi1) != len(umi2):
            raise Exception("various UMI seq length: %s %s" %( umi1, umi2))
            return False
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
                base[0], qual, r1start, r1len, XA[0], *_, q10[0], middle[0], cigartuples[0], reference_start[0], copy_number[0] = reads[0]
                #base[0], qual, r1start, r1len, nmismatch, XA[0], *_, q10[0], middle[0], cigartuples[0], reference_start[0], copy_number[0] = reads[0]
            except TypeError:
                print("Single %s %s " % (reads ,reads[0]))
                continue

            if base[0] in alt:
                if middle[0] is False:
                    nmiddle += 1
                if q10[0] is True:
                    nq10 += 1
                if XA[0]:
                    nmulti += 1
                if copy_number[0] >= 2:
                    nCN_2 += 1
                else:
                    nCN_1 += 1
                noverlap_se += 1
                Alt_set.append(base[0])

        elif len(reads) == 2:
            r1, r2 = reads
            try:
                #base[0], qual, r1start, rlen, nmismatch, XA[0], *_, q10[0], middle[0], cigartuples[0], reference_start[0], copy_number[0] = r1
                #base[1], qual, r1start, rlen, nmismatch, XA[1], *_, q10[1], middle[1], cigartuples[1], reference_start[1], copy_number[1] = r2
                base[0], qual, r1start, rlen, XA[0], *_, q10[0], middle[0], cigartuples[0], reference_start[0], copy_number[0] = r1
                base[1], qual, r1start, rlen, XA[1], *_, q10[1], middle[1], cigartuples[1], reference_start[1], copy_number[1] = r2

            except TypeError:
                print("Double %s %s %s" %(reads, r1, r2))
                continue

            for i, base1 in enumerate(base):
                if base1 in alt:
                    if q10[i] is True:
                        nq10 += 1
                    if middle[i] is False:
                        nmiddle += 1
                    if XA[i]:
                        nmulti += 1
                    if copy_number[i] >= 2:
                        nCN_2 += 1
                    else:
                        nCN_1 += 1
                    noverlap_pe += 0.5
                    
            if base1 in alt:
                Alt_set.append(base1)

        else: # more than 2 reads owing the sam sequence name, aborted!
            nerror += len(reads)

    if umi:
        for name, info in name_dict.items():
            umi_judge = re.search(r':UMI_(\w+)_(\w+)', name)
            if umi_judge:
                if info[0][0] == alt:
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

        numi = 0                            
        for umi1 in UMI_set:
            if umi1 != '' :
                numi += 1

    return nq10, nmiddle, nmulti, int(noverlap_pe), noverlap_se, numi, int(nCN_1), int(nCN_2)