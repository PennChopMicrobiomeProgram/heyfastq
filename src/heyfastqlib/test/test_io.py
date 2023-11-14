from heyfastqlib.io import *

def test_parse_fastq():
    f = [
        "@ab c", "GGCA", "+", "==;G",
        "@d:e:f", "CCGT", "+", "1,4E",
    ]
    seqs = parse_fastq(f)
    assert list(seqs) == [
        ("ab c", "GGCA", "==;G"),
        ("d:e:f", "CCGT", "1,4E"),
    ]

def test_parse_fastq_paired():
    fq1 = ["@a", "TA", "+", "GG", "@b", "CG", "+", "AB"]
    fq2 = ["@a", "AG", "+", "FF", "@b", "TC", "+", "BC"]
    recs = parse_fastq_paired((fq1, fq2))
    assert list(recs) == [
        (("a", "TA", "GG"), ("a", "AG", "FF")),
        (("b", "CG", "AB"), ("b", "TC", "BC")),
    ]

class MockFile:
    def __init__(self):
        self.contents = ""

    def write(self, x):
        self.contents += x

def test_write_fastq():
    f = MockFile()
    seqs = [("a", "CGT", "BBC"), ("b", "TAC", "CCD")]
    write_fastq(f, seqs)
    assert f.contents == "@a\nCGT\n+\nBBC\n@b\nTAC\n+\nCCD\n"

def test_write_fastq_paired():
    f1 = MockFile()
    f2 = MockFile()
    paired_recs = [
        (("a", "CGT", "BBC"), ("a", "ACG", "CCD")),
        (("b", "GTA", "AAB"), ("b", "TAC", "EEF")),
    ]
    write_fastq_paired((f1, f2), paired_recs)
    assert f1.contents == "@a\nCGT\n+\nBBC\n@b\nGTA\n+\nAAB\n"
    assert f2.contents == "@a\nACG\n+\nCCD\n@b\nTAC\n+\nEEF\n"
