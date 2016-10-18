from MrBam.aggregate import aggregate_reads
from argparse import Namespace

def test_aggregate_reads_1():
    "it should aggregate pairs"

    o = Namespace(verbos=False, qual=20)

    reads = (
        ("r1", 'A', 60, 2, 11, 4, 11,  False, True),
        ("r1", 'A', 60, 4, 13, 2, -11, False, True),
        ("r2", 'C', 60, 2, 11, 4, 11,  False, True),
        ("r2", 'C', 60, 4, 13, 2, -11, False, True)
    )

    unique_pairs, *_ = aggregate_reads(o, reads)

    assert len(unique_pairs) == 1

def test_aggregate_reads_2():
    "it should aggregate singles"

    o = Namespace(verbos=False, qual=20)

    reads = (
        ("r1", 'A', 60, 2, 11, 0, 0, False, False),
        ("r2", 'C', 60, 2, 11, 0, 0, False, False),
        ("r3", 'C', 60, 2, 11, 0, 0, True,  False)
    )

    _, unique_single, *_ = aggregate_reads(o, reads)

    assert len(unique_single) == 2

def test_aggregate_reads_3():
    "it should ignore when 3+ reads share the same name"

    o = Namespace(verbos=False, qual=20)

    reads = (
        ("r1", 'A', 60, 2, 11, 2, 9,  False, True),
        ("r1", 'C', 60, 2, 11, 2, -9, False, True),
        ("r1", 'C', 60, 2, 11, 2, 9,  True,  True),
        ("r2", 'C', 60, 2, 11, 0, 0,  True,  False)
    )

    unique_pairs, unique_single, _, nerror, *_ = aggregate_reads(o, reads)

    assert len(unique_pairs) == 0
    assert len(unique_single) == 1
    assert nerror == 3

def test_aggregate_reads_4():
    "it should ignore when base in overlap area inconsistent between two reads"

    o = Namespace(verbos=False, qual=20)

    reads = (
        ("r1", 'A', 60, 2, 11, 4, 11,  False, True),
        ("r1", 'C', 60, 4, 13, 2, -11, False, True),
        ("r2", 'C', 60, 3, 12, 5, 11,  False, True),
        ("r2", 'C', 60, 5, 14, 3, -11, False, True)
    )

    unique_pairs, unique_single, *_, ninconsis = aggregate_reads(o, reads)

    assert len(unique_pairs) == 1
    assert ninconsis == 2
