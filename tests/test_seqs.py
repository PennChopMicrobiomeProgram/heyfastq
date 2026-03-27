from heyfastqlib.seqs import *


def test_kmers():
    seq = "ATGCGCT"
    assert list(kmers(seq, k=4)) == ["ATGC", "TGCG", "GCGC", "CGCT"]
    assert list(kmers(seq, k=5)) == ["ATGCG", "TGCGC", "GCGCT"]


def test_kscore():
    assert kscore("AAAAA", k=4) == 1 / 5
    assert kscore("AAAATAAAAT", k=4) == 5 / 10


def test_kscore_with_ambiguous_bases_matches_original_behavior():
    seq = "ACNGT"
    expected = len(set(kmers(seq, k=2))) / len(seq)
    assert kscore(seq, k=2) == expected


def test_kscore_returns_zero_for_short_sequences():
    assert kscore("AT", k=4) == 0.0
