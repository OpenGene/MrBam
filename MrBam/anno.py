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
    with open(o.query) as fin:
        # vcf file
        for line in fin:
            line = line.rstrip().split()
            chr = line[0]
            pos =int(line[1]) - 1
            try_append(mutation, chr, pos)

            tmp = ''.join((chr,str(pos)))
            variation_info[tmp] = line

    with open(o.output, 'w') as fout:
        for k, sam in enumerate((o.gdna, o.cfdna)):

            variation = deepcopy(mutation) # copy from vcf file
            reads = {} # pos-> read from bam file

            for chr in mutation.keys():
                for pos in mutation[chr]:
                    reads[pos] = []

            LastChr = ''
            LastPos = 0
            for read in sam:
                chr = read.reference_name
                
                if chr != LastChr:
                    LastPos = read.reference_start
                    LastChr = chr
                else:
                    if read.reference_start < LastPos:
                        raise Exception("sam ", sam, " is unsorted!\n")
                    LastPos = read.reference_start

                if chr in variation:
                        pos = variation[chr][0]
                        
                        if read.reference_start <= pos < read.reference_end:
                        # read aligned with first variation site
                            reads[pos].append(read)

                            try:
                                for j in range(1, len(variation[chr])):
                                    if variation[chr][j] < read.reference_end:
                                        reads[variation[chr][j]].append(read)
                                    else:
                                        break
                            except:
                                raise Exception("error in try1\n")

                        elif read.reference_end <= pos:
                            # -- |
                            continue

                        else:
                            # | --
                            variation_info= output(chr, pos, variation_info, o, fout, reads[pos], k, rep_ref)
                            if len(variation[chr]) == 1:
                                del variation[chr]
                                continue
                            else:
                                del variation[chr][0]
                            reads[pos] = []

                            # Reads should be shorter than 200 bases
                            try:
                                for i in range(0,len(variation[chr])) :
                                    # removing variations not aligned with read
                                    pos = variation[chr][0]
                                    if read.reference_start <= pos < read.reference_end:
                                        reads[pos].append(read) 
                                        # the second variation may aligned with vcf file

                                        for j in range(1, len(variation[chr])):
                                            if variation[chr][j] < read.reference_end:
                                                reads[variation[chr][j]].append(read)
                                            else:
                                                break
                                        break

                                    elif read.reference_end <= pos:
                                        # -- |
                                        break

                                    else: # | ---
                                        variation_info= output(chr, pos, variation_info, o, fout, reads[pos], k, rep_ref)
                                        reads[pos] = []
                                        if len(variation[chr]) == 1:
                                            del variation[chr]
                                            break
                                        else:
                                            del variation[chr][0]
                            except:
                                raise Exception("variation ",variation)

            # The last read in one chr aligned with mutation site, leaving mutation info undealed
            for chr in variation:
                for pos in variation[chr]:
                    if reads[pos] != []:
                        variation_info= output(chr, pos, variation_info, o, fout, reads[pos], k, rep_ref)
                    else:
                        raise Exception("reads issue ",chr, pos, variation)


def output(chr, pos, variation_info, o, fout, reads_pos, k, rep_ref):
    tmp = ''.join((chr,str(pos)))
    line = variation_info[tmp]
    ref = line[3]
    alt = line[4]

    new_reads = get_infor(o, reads_pos, pos, ref)
    name_dict, unique_pairs, unique_single= snv_mut(new_reads, o, None if o.fast else pad_softclip(sam))
    mor, mnr, msr, oor, onr, osr, moa, mna, msa, ooa, ona, osa, _ = count_different_type(o, unique_pairs, unique_single, alt, ref)

    if k == 0:
        line[-3] += ':UDP'

    line[-2] += ':' + ','.join(map(str, (mor, mnr, msr, oor, onr, osr, moa, mna, msa, ooa, ona, osa)))

    (nq10_a, nterminal_a, nmulti_a, noverlap_pe_a, noverlap_se_a, numi_a, nCN_1_a, nCN_2_a, ave_mapqual_a ) = extra_info(o, name_dict, ref, alt)
    
    repeat = repeat_area(rep_ref, chr, pos + 1, o)
    extra_2 = ','.join(map(str, (repeat, nq10_a, nterminal_a, nmulti_a, noverlap_pe_a, noverlap_se_a, numi_a, nCN_1_a, nCN_2_a, ave_mapqual_a)))
    line.append(extra_2)

    variation_info[tmp] = line
    if k == 1:
        print(file=fout, sep='\t', *line)
        fout.flush()
        del variation_info[tmp]
    
    return variation_info


def repeat_area(rep_ref, chr, pos, o):
    if o.repeat is None:
        return 'NA'

    if chr in rep_ref:
        if int(rep_ref[chr][0] <= pos < int(rep_ref[chr][1]) + 1):
            return True
        else:
            return False                             
    return False