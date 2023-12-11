from heyfastqlib.read import *


def test_trim():
    assert trim(Read("myseq", "ACGTAC", "123456"), 4) == Read("myseq", "ACGT", "1234")


def test_length():
    assert length(Read("a", "ATCGC", "12345")) == 5


def test_kscore_ok():
    obs = Read("a", "AAAAC", "12345")  # kscore = 2 / 5 = 0.2
    assert not kscore_ok(obs, min_kscore=0.5)
    assert kscore_ok(obs, min_kscore=0.1)


def test_length_ok():
    obs = Read("g", "ACTTACT", "1234567")
    assert length_ok(obs, 5)
    assert not length_ok(obs, 10)
