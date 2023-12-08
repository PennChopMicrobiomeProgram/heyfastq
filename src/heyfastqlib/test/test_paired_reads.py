from heyfastqlib.paired_reads import *

def test_trim_fixed():
    assert trim_fixed(("myseq", "ACGTAC", "123456"), 4) == \
        ("myseq", "ACGT", "1234")

def test_kscore_ok():
    obs = ("a", "AAAAC", "12345") # kscore = 2 / 5 = 0.2
    assert not kscore_ok(obs, min_kscore=0.5)
    assert kscore_ok(obs, min_kscore=0.1)

reads = [
    (("a", "AAAG", "1234"), ("a", "CTTT", "1234")),
    (("b", "ACGTACGT", "12345678"), ("b", "ACG", "123")),
]

def test_map_paired():
    assert list(map_paired(reads, read_length)) == [(4, 4), (8, 3)]

def test_filter_paired():
    assert list(filter_paired(reads, compare_length, threshold = 4)) == \
        reads[:1]

    assert list(
        filter_paired(reads, compare_length, requirement=any, threshold = 4)
    ) == reads
