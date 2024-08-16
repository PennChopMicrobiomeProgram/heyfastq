from heyfastqlib.io import *


def test_parse_fastq():
    f = [
        "@ab c",
        "GGCA",
        "+",
        "==;G",
        "@d:e:f",
        "CCGT",
        "+",
        "1,4E",
    ]
    reads = parse_fastq(f)
    assert list(reads) == [
        Read("ab c", "GGCA", "==;G"),
        Read("d:e:f", "CCGT", "1,4E"),
    ]


def test_parse_fastq_paired():
    fq1 = ["@a", "TA", "+", "GG", "@b", "CG", "+", "AB"]
    fq2 = ["@a", "AG", "+", "FF", "@b", "TC", "+", "BC"]
    recs = parse_fastq_paired((fq1, fq2))
    assert list(recs) == [
        (Read("a", "TA", "GG"), Read("a", "AG", "FF")),
        (Read("b", "CG", "AB"), Read("b", "TC", "BC")),
    ]


class MockFile:
    def __init__(self):
        self.contents = ""

    def write(self, x):
        self.contents += x


def test_write_fastq():
    f = MockFile()
    reads = [Read("a", "CGT", "BBC"), Read("b", "TAC", "CCD")]
    write_fastq(f, reads)
    assert f.contents == "@a\nCGT\n+\nBBC\n@b\nTAC\n+\nCCD\n"


def test_write_fastq_paired():
    f1 = MockFile()
    f2 = MockFile()
    paired_recs = [
        (Read("a", "CGT", "BBC"), Read("a", "ACG", "CCD")),
        (Read("b", "GTA", "AAB"), Read("b", "TAC", "EEF")),
    ]
    write_fastq_paired((f1, f2), paired_recs)
    assert f1.contents == "@a\nCGT\n+\nBBC\n@b\nGTA\n+\nAAB\n"
    assert f2.contents == "@a\nACG\n+\nCCD\n@b\nTAC\n+\nEEF\n"


def test_parse_seq_ids():
    f = ["Id1\n", "\tId2|345 678  \n", "   ", "   # a comment", "  id3"]
    assert list(parse_seq_ids(f)) == ["Id1", "Id2|345", "id3"]
