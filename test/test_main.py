from MrBam.main import init
from helper import make_bam
from argparse import Namespace
from pysam import AlignmentFile
from pytest import raises

def test_init_1():
    "it should raise exception if neither of cfdna and gdna are not provided"
    o = Namespace(query="test.vcf", cfdna=None, gdna=None, output=None)
    with raises(Exception):
        init(o)

def test_init_2(tmpdir):
    "it should open sam file if provided"

    make_bam(tmpdir.strpath, """
             123456789_123456789_
        r1 + ...........
        r1 -       ......*....
        r2 +    .........*.
        r2 -        .....*.......
    """)

    o = Namespace(query="test.vcf", cfdna=tmpdir.join("test.bam").strpath, gdna=None, output=None)

    init(o)

    assert isinstance(o.cfdna, AlignmentFile)
    assert o.gdna == None

def test_init_3(tmpdir):
    "it should generate a proper output file name if not provided"

    make_bam(tmpdir.strpath, """
             123456789_123456789_
        r1 + ...........
        r1 -       ......*....
        r2 +    .........*.
        r2 -        .....*.......
    """)

    o = Namespace(query="test.vcf", cfdna=tmpdir.join("test.bam").strpath, gdna=None, output=None)

    init(o)

    assert o.output != None
