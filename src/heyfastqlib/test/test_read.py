from heyfastqlib.read import (
    Read,
    length,
    qvals,
    trim,
    kscore_ok,
    length_ok,
    trim_moving_average,
    trim_ends,
    seq_id_ok,
)


def test_length():
    assert length(Read("a", "ATCGC", "12345")) == 5


def test_qvals():
    assert qvals(Read("a", "ACGT", "!.:@F")) == [0, 13, 25, 31, 37]


def test_trim():
    a = Read("myseq", "ACGTAC", "123456")
    assert trim(a, end_idx=4) == Read("myseq", "ACGT", "1234")


def test_trim_returns_copy():
    a = Read("a", "ACGG", "FFFF")
    assert length(trim(a, 0, 2)) == 2
    assert length(a) == 4


def test_kscore_ok():
    obs = Read("a", "AAAAC", "12345")  # kscore = 2 / 5 = 0.2
    assert not kscore_ok(obs, min_kscore=0.5)
    assert kscore_ok(obs, min_kscore=0.1)


def test_length_ok():
    obs = Read("g", "ACTTACT", "1234567")
    assert length_ok(obs, 5)
    assert not length_ok(obs, 10)


def test_trim_moving_average():
    r = Read("a", "ACGTACGTAAAAAA", "FFFFFFFF......")
    assert trim_moving_average(r, 3, 15) == Read("a", "ACGTACGT", "FFFFFFFF")


def test_trim_moving_average_endcaps():
    # idx   0  1  2  3 4 5
    # qual 40 40 40 40 0 0
    # ave  40 30 20 10 0 0
    # w25   +  +  -  - - -
    # q25   +  +  +  + - -
    # Window falls below Q25 at idx 2
    # Within window, keep idx 2, 3
    # Trim at 4
    a = Read("a", "ACGTAAAAAAAA", "IIII!!!!!!!!")
    assert trim_moving_average(a, 4, 25) == trim(a, 0, 4)

    # idx   0  1  2  3  4  5 6 7 8 9
    # qual 40 40 40 40  0 40 0 0 0 0
    # ave  40 30 30 20 10 10 0
    # w25   +  +  +  -  -  - -
    # q25   +  +  +  +  -  + - - - -
    # Window falls below Q25 at idx 3
    # Within window, keep 3, 4, 5
    # Trim at 6
    b = Read("b", "CGTTCCCCCCCC", "IIII!I!!!!!!")
    assert trim_moving_average(b, 4, 25) == trim(b, 0, 6)

    # idx   0  1  2  3  4  5  6 7 8 9
    # qual 40 40 40 40  0  0 40 0 0 0
    # ave  40 30 20 20 10 10 10
    # w25   +  +  -  -  -  -  -
    # q25   +  +  +  +  -  -  +
    # Window falls below Q25 at idx 2
    # Within window, keep 2, 3
    # Trim at 4
    c = Read("c", "GCGGACGTCGGG", "IIII!!I!!!!!")
    assert trim_moving_average(c, 4, 25) == trim(c, 0, 4)

    # idx   0  1  2  3  4  5  6 7 8 9
    # qual 40 40 40 40  0  0 40 0 0 0
    # ave  40 30 20 20 10 10 10 0
    # w15   +  +  +  +  -  -  - -
    # q15   +  +  +  +  -  -  + -
    # Window falls below Q15 at idx 4
    # Within window, keep 4, 5, 6
    # Trim at 7
    d = Read("d", "GCGGACGTCGGG", "IIII!!I!!!!!")
    assert trim_moving_average(d, 4, 15) == trim(d, 0, 7)


def test_trim_ends():
    a = Read("bf", "ACGTACGT", "&!FFF!.&")
    # idx  0 1  2  3  4 5  6 7
    # qual 5 0 37 37 37 0 13 5
    # q20  - -  +  +  + -  - - Trim 2:5
    # q6   - -  +  +  + -  + - Trim 2:7
    assert trim_ends(a, 20, 20) == trim(a, 2, 5)
    assert trim_ends(a, 6, 6) == trim(a, 2, 7)


def test_seq_id_ok():
    ids = ["cd", "b"]
    assert seq_id_ok(Read("b", "GGC", "FFF"), ids)
    assert not seq_id_ok(Read("a", "TCG", "FFF"), ids)
    assert not seq_id_ok(Read("cd", "TCG", "!!!"), ids, remove=True)
