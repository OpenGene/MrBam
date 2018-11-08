
def explain():
    print("example result\t","1026,202,0,344,338,189,0,682,500,60,194")
    print("First Column\t", "the number of reads support mutation")
    print("Second Column\t", "the number of reads owing mutation located at the 5' or 3' of the read, distance cutoff: 10%")
    print("Third Column\t", "the number of reads owing mutation and are multiple alignments by BWA")
    print("Fourth Column\t", "the number of molecules owing mutation located the overlap area of pair-end reads")
    print("Fifth Column\t", "the number of molecules that are single end reads or one of read owing mutation while the other not")
    print("Sixth Column\t", "the number of deduplicated UMI")
    print("Seventh Column\t", "the number of molecules that got copy number equal to 1")
    print("Eighth Column\t", "the number of molecules that got copy number more than 1")
    print('Ninth Column\t', "the number of molecules without XA, mapping quality more than 0 and mutation site located at the middle of read")
    print("Tenth Column\t", "the average map quality of reads support mutation")
    print("Eleventh Column\t", "the average insert size of reads support mutation")
    exit()

    