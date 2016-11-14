from MrBam.anno import anno
from argparse import Namespace
from helper import make_bam
from pysam import AlignmentFile

def test_anno_1(tmpdir):
    "test --simple version"

    make_bam(tmpdir.strpath, """
             123456789_123456789_12
        r1 + ...........
        r1 -      ......*....
        r2 +   .........*.
        r2 -       .....*.......
        r3 +       ...........
        r3 -            ....*......
        r4 +       ...........
        r4 -            ...........
             123456789_123456789_12
    """)

    sam = AlignmentFile(tmpdir.join("test.bam").strpath)

    o = Namespace(
        query  = tmpdir.join("test.vcf").strpath,
        cfdna  = sam,
        gdna   = None,
        simple = True,
        verbos = False,
        fast   = False,
        qual   = 20,
        output = tmpdir.join("test_MrBam.vcf").strpath,
        allow_inconsist = False,
        mismatch_limit  = -1,
    )

    anno(o)

    for i in open(tmpdir.join("test_MrBam.vcf").strpath):
        if i.startswith('#'):
            continue

        i = i.split('\t')

        if i[1] == '12':
            a,b,c,d = i[-1].split(':')[-1].strip().split(',')
            assert a == '0'
            assert b == '0'
            assert c == '1'
            assert d == '1'
        elif i[1] == '16':
            a,b,c,d = i[-1].split(':')[-1].strip().split(',')
            assert a == '0'
            assert b == '0'
            assert c == '0'
            assert d == '0'
        else:
            raise Exception("unexpected variant call")
