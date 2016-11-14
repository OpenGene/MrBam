from MrBam.bam import get_reads, pad_softclip
from helper import make_bam
from argparse import Namespace
from pysam import AlignmentFile

def test_get_reads_1(tmpdir):
    "it should get all but only the reads that covers the given position"

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

    o = Namespace(verbos=False, mismatch_limit=-1)
    sam = AlignmentFile(tmpdir.join("test.bam").strpath)

    assert sum( 1 for _ in get_reads(o, sam, 'ref', '4') )  == 2
    assert sum( 1 for _ in get_reads(o, sam, 'ref', '12') ) == 7
    assert sum( 1 for _ in get_reads(o, sam, 'ref', '20') ) == 2

def test_get_reads_2(tmpdir):
    "it should read properties correctly"

    make_bam(tmpdir.strpath, """
             123456789_123
        r1 + ...*.......
        r1 -   .*.........
    """)

    o = Namespace(verbos=False, mismatch_limit=-1)
    sam = AlignmentFile(tmpdir.join("test.bam").strpath)

    r = next(get_reads(o, sam, 'ref', '4'))

    assert r[0] == "r1" # name
    assert r[3] == 0 # 0-based pos
    assert r[4] == 11 # length
    assert r[5] == -1 # mismatch, not caculated
    assert r[6] == 2 # mate pos
    assert r[7] == 13 # template length
    assert r[8] == False # is_reverse
    assert r[9] == True # paired and mapped

def test_pad_softclip_1(tmpdir):
    "it should memorize the result"

    make_bam(tmpdir.strpath, """
        r1 + __.*.......
        r1 -   .*.......__
    """)

    o = Namespace(verbos=False, mismatch_limit=-1)
    sam = AlignmentFile(tmpdir.join("test.bam").strpath)

    a = pad_softclip(sam)
    b = pad_softclip(sam)

    assert a is b

def test_pad_softclip_2(tmpdir):
    "it should ignore more than two reads which share the same name"

    make_bam(tmpdir.strpath, """
        r1 + __.*.......
        r1 -   .*.......__
        r1 -   .*.......__
        r2 +   .*.......__
        r2 -   .*.......__
    """)

    o = Namespace(verbos=False, mismatch_limit=-1)
    sam = AlignmentFile(tmpdir.join("test.bam").strpath)

    adjusted_pos = pad_softclip(sam)

    assert sum(1 for startpos, length in adjusted_pos.values() if startpos != -1) == 1

def test_pad_softclip_3(tmpdir):
    "it should pad softclipped bases"

    make_bam(tmpdir.strpath, """
             123456789_123
        r1 + __.*.......
        r1 -   .*.........
        r2 - ...*.......
        r2 +   .*.......__
    """)

    o = Namespace(verbos=False, mismatch_limit=-1)
    sam = AlignmentFile(tmpdir.join("test.bam").strpath)

    adjusted_pos = pad_softclip(sam)

    assert adjusted_pos["r1"] == (0, 13) # 0-based position
    assert adjusted_pos["r2"] == (0, 13)
