from heyfastqlib.read import Read
from heyfastqlib.diff import FastqDiff, diff_offset

def test_diff_same():
    a1 = "ATCGCGTCA"
    a2 = "ATCGCGTCA"
    assert diff_offset(a1, a2) == 0

def test_diff_mismatch():
    a1 = "ATCCCGTCA"
    a2 = "ATCGCGTCA"
    assert diff_offset(a1, a2) == 0

def test_diff_trim_right():
    a1 = "ATCGCGTCA"
    a2 = "ATCGCG"
    assert diff_offset(a1, a2) == 0
    
def test_diff_trim_left():
    a1 = "ATCGCGTCA"
    a2 =   "CGCGTCA"
    assert diff_offset(a1, a2) == 2

def test_diff_trim_both():
    a1 = "ATCGCGTCA"
    a2 =  "TCGCG"
    assert diff_offset(a1, a2) == 1

def test_diff_extend_right():
    a1 = "ATCGCG"
    a2 = "ATCGCGTCA"
    assert diff_offset(a1, a2) == 0

def test_diff_extend_left():
    a1 =  "TCGCGTCA"
    a2 = "ATCGCGTCA"
    assert diff_offset(a1, a2) == -1

def test_diff_extend_both():
    a1 =   "CGCGT"
    a2 = "ATCGCGTCA"
    assert diff_offset(a1, a2) == -2

def test_diff_format():
    d = FastqDiff("read1", "GCTCA", "AGCCCA", -1)
    assert list(d.format()) == [
        "read1",
        "   GCTCA",
        "   || ||",
        "  AGCCCA",
    ]

def test_diff_from_reads():
    a = Read("read1", "CGCGT", ".....")
    b = Read("read1", "ATCGCGTCA", ".........")
    assert FastqDiff.from_reads(a, b) == \
        FastqDiff("read1", "CGCGT", "ATCGCGTCA", -2)
