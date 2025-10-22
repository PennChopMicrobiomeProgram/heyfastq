from io import StringIO

from heyfastqlib.io import parse_fastq, parse_fastq_single, parse_seq_ids, write_fastq
from heyfastqlib.read import Read


def make_fastq(lines):
    return StringIO("\n".join(lines) + "\n")


def test_parse_fastq():
    handle = make_fastq(
        [
            "@ab c",
            "GGCA",
            "+",
            "==;G",
            "@d:e:f",
            "CCGT",
            "+",
            "1,4E",
        ]
    )
    reads = parse_fastq_single(handle)
    assert list(reads) == [
        Read("ab c", "GGCA", "==;G"),
        Read("d:e:f", "CCGT", "1,4E"),
    ]


def test_parse_fastq_paired():
    fq1 = make_fastq(["@a", "TA", "+", "GG", "@b", "CG", "+", "AB"])
    fq2 = make_fastq(["@a", "AG", "+", "FF", "@b", "TC", "+", "BC"])
    recs = parse_fastq((fq1, fq2))
    assert list(recs) == [
        (Read("a", "TA", "GG"), Read("a", "AG", "FF")),
        (Read("b", "CG", "AB"), Read("b", "TC", "BC")),
    ]


def test_write_fastq():
    dest = StringIO()
    reads = [Read("a", "CGT", "BBC"), Read("b", "TAC", "CCD")]
    write_fastq((dest,), iter(reads))
    assert dest.getvalue() == "@a\nCGT\n+\nBBC\n@b\nTAC\n+\nCCD\n"


def test_write_fastq_paired():
    dest1 = StringIO()
    dest2 = StringIO()
    paired_recs = [
        (Read("a", "CGT", "BBC"), Read("a", "ACG", "CCD")),
        (Read("b", "GTA", "AAB"), Read("b", "TAC", "EEF")),
    ]
    write_fastq((dest1, dest2), iter(paired_recs))
    assert dest1.getvalue() == "@a\nCGT\n+\nBBC\n@b\nGTA\n+\nAAB\n"
    assert dest2.getvalue() == "@a\nACG\n+\nCCD\n@b\nTAC\n+\nEEF\n"


def test_parse_seq_ids():
    handle = StringIO("Id1\n\tId2|345 678  \n   \n   # a comment\n  id3")
    assert list(parse_seq_ids(handle)) == ["Id1", "Id2|345", "id3"]
