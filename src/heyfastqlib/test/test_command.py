import gzip
from heyfastqlib.command import *


def test_command(tmp_path):
    in1 = tmp_path / "input_1.fastq"
    with open(in1, "w") as f:
        f.write("@a\nCGTT\n+\n;=GG\n@b\nACTG\n+\nGGGG\n")
    in2 = tmp_path / "input_2.fastq"
    with open(in2, "w") as f:
        f.write("@a\nAACG\n+\n;=GG\n@b\nCAGT\n+\nGGGG\n")
    out1 = tmp_path / "output_1.fastq"
    out2 = tmp_path / "output_2.fastq"
    heyfastq_main(
        [
            "trim-fixed",
            "--length",
            "2",
            "--input",
            str(in1),
            str(in2),
            "--output",
            str(out1),
            str(out2),
        ]
    )
    with open(out1) as f:
        assert f.read() == "@a\nCG\n+\n;=\n@b\nAC\n+\nGG\n"
    with open(out2) as f:
        assert f.read() == "@a\nAA\n+\n;=\n@b\nCA\n+\nGG\n"


def test_gzip_command(tmp_path):
    in1 = tmp_path / "input_1.fastq.gz"
    with gzip.open(in1, "wt") as f:
        f.write("@a\nCGTT\n+\n;=GG\n@b\nACTG\n+\nGGGG\n")
    in2 = tmp_path / "input_2.fastq.gz"
    with gzip.open(in2, "wt") as f:
        f.write("@a\nAACG\n+\n;=GG\n@b\nCAGT\n+\nGGGG\n")
    out1 = tmp_path / "output_1.fastq.gz"
    out2 = tmp_path / "output_2.fastq.gz"
    heyfastq_main(
        [
            "trim-fixed",
            "--length",
            "2",
            "--input",
            str(in1),
            str(in2),
            "--output",
            str(out1),
            str(out2),
        ]
    )
    with gzip.open(out1, "rt") as f:
        assert f.read() == "@a\nCG\n+\n;=\n@b\nAC\n+\nGG\n"
    with gzip.open(out2, "rt") as f:
        assert f.read() == "@a\nAA\n+\n;=\n@b\nCA\n+\nGG\n"


def test_trim_qual_command(tmp_path):
    in1 = tmp_path / "input_1.fastq"
    with open(in1, "w") as f:
        f.write("@a\nAACGTACGTGGGGGGGG\n+\n&55555555&&&&&&&&")
    in2 = tmp_path / "input_2.fastq"
    with open(in2, "w") as f:
        f.write("@a\nCGTTCGTTAAAAAAAA\n+\n55555555!!!!!!!!")
    out1 = tmp_path / "output_1.fastq"
    out2 = tmp_path / "output_2.fastq"
    heyfastq_main(
        [
            "trim-qual",
            "--window-width",
            "4",
            "--window-threshold",
            "7",
            "--start-threshold",
            "6",
            "--min-length",
            "4",
            "--input",
            str(in1),
            str(in2),
            "--output",
            str(out1),
            str(out2),
        ]
    )
    with open(out1) as f:
        # Final Q20 base passes
        # mean([20, 5, 5, 5]) = 8.75
        assert f.read() == "@a\nACGTACGT\n+\n55555555\n"
    with open(out2) as f:
        # Final Q20 base fails, but is added in extension
        # mean([20, 0, 0, 0]) = 5
        assert f.read() == "@a\nCGTTCGTT\n+\n55555555\n"


in1_kscore = """\
@a
AAAATAAAAAAAAAA
+
===============
@b
GCTAGCTAGCATGCATCTA
+
===================
"""

in2_kscore = """\
@a
GCTACGATCAGTACGAT
+
=================
@b
GCTGAGCTACGGTC
+
==============
"""


def test_filter_kscore_command(tmp_path):
    in1 = tmp_path / "input_1.fastq"
    with open(in1, "w") as f:
        f.write(in1_kscore)
    in2 = tmp_path / "input_2.fastq"
    with open(in2, "w") as f:
        f.write(in2_kscore)
    out1 = tmp_path / "output_1.fastq"
    out2 = tmp_path / "output_2.fastq"
    heyfastq_main(
        [
            "filter-kscore",
            "--min-kscore",
            "0.55",
            "--input",
            str(in1),
            str(in2),
            "--output",
            str(out1),
            str(out2),
        ]
    )
    with open(out1) as f:
        assert f.read() == "@b\nGCTAGCTAGCATGCATCTA\n+\n===================\n"
    with open(out2) as f:
        assert f.read() == "@b\nGCTGAGCTACGGTC\n+\n==============\n"


in1_length = """\
@a
ACGTACGTACGT
+
123456789012
@b
AAGGC
+
12345
"""

in2_length = """\
@a
AGGTCGTCTAAC
+
123456789012
@b
AGCTGCTACGCTA
+
1234567890123
"""


def test_filter_length_command(tmp_path):
    in1 = tmp_path / "input_1.fastq"
    with open(in1, "w") as f:
        f.write(in1_length)
    in2 = tmp_path / "input_2.fastq"
    with open(in2, "w") as f:
        f.write(in2_length)
    out1 = tmp_path / "output_1.fastq"
    out2 = tmp_path / "output_2.fastq"
    heyfastq_main(
        [
            "filter-length",
            "--length",
            "6",
            "--input",
            str(in1),
            str(in2),
            "--output",
            str(out1),
            str(out2),
        ]
    )
    with open(out1) as f:
        assert f.read() == "@a\nACGTACGTACGT\n+\n123456789012\n"
    with open(out2) as f:
        assert f.read() == "@a\nAGGTCGTCTAAC\n+\n123456789012\n"
