from heyfastqlib.seqs import *


def test_kmers():
    seq = "ATGCGCT"
    assert list(kmers(seq, k=4)) == ["ATGC", "TGCG", "GCGC", "CGCT"]
    assert list(kmers(seq, k=5)) == ["ATGCG", "TGCGC", "GCGCT"]


def test_kscore():
    assert kscore("AAAAA", k=4) == 1 / 5
    assert kscore("AAAATAAAAT", k=4) == 5 / 10
