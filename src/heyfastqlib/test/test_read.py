from heyfastqlib.read import (
    Read,
    length,
    qvals,
    trim,
    kscore_ok,
    length_ok,
    trim_moving_average,
)


def test_length():
    assert length(Read("a", "ATCGC", "12345")) == 5


def test_qvals():
    assert qvals(Read("a", "ACGT", "!.:@F")) == [0, 13, 25, 31, 37]


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
    assert trim_moving_average(a, 4, 25) == Read("a", "ACGT", "IIII")

    # idx   0  1  2  3  4  5 6 7 8 9
    # qual 40 40 40 40  0 40 0 0 0 0
    # ave  40 30 30 20 10 10 0
    # w25   +  +  +  -  -  - -
    # q25   +  +  +  +  -  + - - - -
    # Window falls below Q25 at idx 3
    # Within window, keep 3, 4, 5
    # Trim at 6
    b = Read("b", "CGTTCCCCCCCC", "IIII!I!!!!!!")
    assert trim_moving_average(b, 4, 25) == Read("b", "CGTTCC", "IIII!I")

    # idx   0  1  2  3  4  5  6 7 8 9
    # qual 40 40 40 40  0  0 40 0 0 0
    # ave  40 30 20 20 10 10 10
    # w25   +  +  -  -  -  -  -
    # q25   +  +  +  +  -  -  +
    # Window falls below Q25 at idx 2
    # Within window, keep 2, 3
    # Trim at 4
    c = Read("c", "GCGGACGTCGGG", "IIII!!I!!!!!")
    assert trim_moving_average(c, 4, 25) == Read("c", "GCGG", "IIII")

    # idx   0  1  2  3  4  5  6 7 8 9
    # qual 40 40 40 40  0  0 40 0 0 0
    # ave  40 30 20 20 10 10 10 0
    # w15   +  +  +  +  -  -  - -
    # q15   +  +  +  +  -  -  + -
    # Window falls below Q15 at idx 4
    # Within window, keep 4, 5, 6
    # Trim at 7
    d = Read("d", "GCGGACGTCGGG", "IIII!!I!!!!!")
    assert trim_moving_average(d, 4, 15) == Read("d", "GCGGACG", "IIII!!I")
