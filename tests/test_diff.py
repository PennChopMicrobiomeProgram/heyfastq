from heyfastqlib.read import Read
from heyfastqlib.diff import ReadAligned, diff_offset


def test_diff_same():
    a1 = "ATCGCGTCA"
    a2 = "ATCGCGTCA"
    assert diff_offset(a1, a2) == 0
    d = ReadAligned("x", a1, a2, 0)
    assert d.left_trim == ""
    assert d.aligned1 == a1
    assert d.right_trim == ""
    assert d.left_extend == ""
    assert d.aligned2 == a2
    assert d.right_extend == ""


def test_diff_mismatch():
    a1 = "ATCCCGTCA"
    a2 = "ATCGCGTCA"
    assert diff_offset(a1, a2) == 0
    d = ReadAligned("x", a1, a2, 0)
    assert d.left_trim == ""
    assert d.aligned1 == a1
    assert d.right_trim == ""
    assert d.left_extend == ""
    assert d.aligned2 == a2
    assert d.right_extend == ""


def test_diff_trim_right():
    a1 = "ATCGCGTCA"
    a2 = "ATCGCG"
    assert diff_offset(a1, a2) == 0
    d = ReadAligned("x", a1, a2, 0)
    assert d.left_trim == ""
    assert d.aligned1 == "ATCGCG"
    assert d.right_trim == "TCA"
    assert d.left_extend == ""
    assert d.aligned2 == a2
    assert d.right_extend == ""


def test_diff_trim_left():
    a1 = "ATCGCGTCA"
    a2 =   "CGCGTCA"  # fmt: skip
    assert diff_offset(a1, a2) == 2
    d = ReadAligned("x", a1, a2, 2)
    assert d.left_trim == "AT"
    assert d.aligned1 == "CGCGTCA"
    assert d.right_trim == ""
    assert d.left_extend == ""
    assert d.aligned2 == a2
    assert d.right_extend == ""


def test_diff_trim_both():
    a1 = "ATCGCGTCA"
    a2 =  "TCGCG"  # fmt: skip
    assert diff_offset(a1, a2) == 1
    d = ReadAligned("x", a1, a2, 1)
    assert d.left_trim == "A"
    assert d.aligned1 == "TCGCG"
    assert d.right_trim == "TCA"
    assert d.left_extend == ""
    assert d.aligned2 == a2
    assert d.right_extend == ""


def test_diff_extend_right():
    a1 = "ATCGCG"
    a2 = "ATCGCGTCA"
    assert diff_offset(a1, a2) == 0
    d = ReadAligned("x", a1, a2, 0)
    assert d.left_trim == ""
    assert d.aligned1 == a1
    assert d.right_trim == ""
    assert d.left_extend == ""
    assert d.aligned2 == "ATCGCG"
    assert d.right_extend == "TCA"


def test_diff_extend_left():
    a1 =  "TCGCGTCA"  # fmt: skip
    a2 = "ATCGCGTCA"
    assert diff_offset(a1, a2) == -1
    d = ReadAligned("x", a1, a2, -1)
    assert d.left_trim == ""
    assert d.aligned1 == a1
    assert d.right_trim == ""
    assert d.left_extend == "A"
    assert d.aligned2 == "TCGCGTCA"
    assert d.right_extend == ""


def test_diff_extend_both():
    a1 =   "CGCGT"  # fmt: skip
    a2 = "ATCGCGTCA"
    assert diff_offset(a1, a2) == -2
    d = ReadAligned("x", a1, a2, -2)
    assert d.left_trim == ""
    assert d.aligned1 == a1
    assert d.right_trim == ""
    assert d.left_extend == "AT"
    assert d.aligned2 == "CGCGT"
    assert d.right_extend == "CA"


def test_diff_empty():
    assert diff_offset("ACGT", "") == 0
    assert diff_offset("", "AGTC") == 0
    d = ReadAligned("x", "ACGT", "", 0)
    assert d.left_trim == ""
    assert d.aligned1 == ""
    assert d.right_trim == "ACGT"
    assert d.left_extend == ""
    assert d.aligned2 == ""
    assert d.right_extend == ""


def test_diff_format_alignment():
    d = ReadAligned("read1", "GCTCA", "AGCCCA", -1)
    assert d.format_alignment() == "read1\n   GCTCA\n   || ||\n  AGCCCA\n"


def test_diff_format_alignmment_empty():
    d = ReadAligned("a", "GCATCG", "", 0)
    assert d.format_alignment() == "a\n  GCATCG\n        \n        \n"


def test_diff_from_reads():
    a = Read("x", "CGCG", ".....")
    b = Read("x", "ATCGCGTC", ".........")
    assert ReadAligned.from_reads(a, b) == ReadAligned("x", "CGCG", "ATCGCGTC", -2)

def test_diff_format_tsv():
    d = ReadAligned("x", "ACGTG", "CAA", 1)
    assert d.format_tsv() == "x\ttrim-both-mismatch\tA\t\tG\t\t2,G,A;3,T,A\n"
