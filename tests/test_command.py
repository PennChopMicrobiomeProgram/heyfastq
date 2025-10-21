import gzip
import shutil
from pathlib import Path

from heyfastqlib.command import fastq_io_parser, heyfastq_main

DATA_DIR = Path(__file__).parent / "data"

if not any(action.dest == "threads" for action in fastq_io_parser._actions):
    fastq_io_parser.add_argument("--threads", type=int, default=1)


def copy_data(tmp_path, filename):
    src = DATA_DIR / filename
    dest = tmp_path / filename
    shutil.copyfile(src, dest)
    return dest


def read_expected(filename):
    return (DATA_DIR / filename).read_text()


def test_trim_fixed_command(tmp_path):
    in1 = copy_data(tmp_path, "trim_fixed_input_1.fastq")
    in2 = copy_data(tmp_path, "trim_fixed_input_2.fastq")
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

    assert out1.read_text() == "@a\nCG\n+\n;=\n@b\nAC\n+\nGG\n"
    assert out2.read_text() == "@a\nAA\n+\n;=\n@b\nCA\n+\nGG\n"


def test_gzip_command(tmp_path):
    content1 = read_expected("trim_fixed_input_1.fastq")
    content2 = read_expected("trim_fixed_input_2.fastq")
    in1 = tmp_path / "input_1.fastq.gz"
    in2 = tmp_path / "input_2.fastq.gz"

    with gzip.open(in1, "wt") as handle:
        handle.write(content1)
    with gzip.open(in2, "wt") as handle:
        handle.write(content2)

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

    with gzip.open(out1, "rt") as handle:
        assert handle.read() == "@a\nCG\n+\n;=\n@b\nAC\n+\nGG\n"
    with gzip.open(out2, "rt") as handle:
        assert handle.read() == "@a\nAA\n+\n;=\n@b\nCA\n+\nGG\n"


def test_trim_qual_command(tmp_path):
    in1 = copy_data(tmp_path, "trim_qual_input_1.fastq")
    in2 = copy_data(tmp_path, "trim_qual_input_2.fastq")
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

    assert out1.read_text() == "@a\nACGTACGT\n+\n55555555\n"
    assert out2.read_text() == "@a\nCGTTCGTT\n+\n55555555\n"


def test_filter_kscore_command(tmp_path):
    in1 = copy_data(tmp_path, "filter_kscore_input_1.fastq")
    in2 = copy_data(tmp_path, "filter_kscore_input_2.fastq")
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

    assert out1.read_text() == "@b\nGCTAGCTAGCATGCATCTA\n+\n===================\n"
    assert out2.read_text() == "@b\nGCTGAGCTACGGTC\n+\n==============\n"


def test_filter_length_command(tmp_path):
    in1 = copy_data(tmp_path, "filter_length_input_1.fastq")
    in2 = copy_data(tmp_path, "filter_length_input_2.fastq")
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

    assert out1.read_text() == "@a\nACGTACGTACGT\n+\n123456789012\n"
    assert out2.read_text() == "@a\nAGGTCGTCTAAC\n+\n123456789012\n"


def test_filter_seq_ids_command(tmp_path):
    seqids = copy_data(tmp_path, "filter_seqids_ids.txt")
    in1 = copy_data(tmp_path, "filter_seqids_input.fastq")
    in2 = copy_data(tmp_path, "filter_seqids_input_2.fastq")
    out1 = tmp_path / "output_1.fastq"
    out2 = tmp_path / "output_2.fastq"

    heyfastq_main(
        [
            "filter-seqids",
            str(seqids),
            "--input",
            str(in1),
            str(in2),
            "--output",
            str(out1),
            str(out2),
        ]
    )

    assert out1.read_text() == "@b\nGTCC\n+\n5678\n"
    assert out2.read_text() == "@b\nCCGG\n+\n####\n"


def test_subsample_command(tmp_path):
    in1 = copy_data(tmp_path, "subsample_input_1.fastq")
    in2 = copy_data(tmp_path, "subsample_input_2.fastq")
    out1 = tmp_path / "output_1.fastq"
    out2 = tmp_path / "output_2.fastq"

    heyfastq_main(
        [
            "subsample",
            "--n",
            "2",
            "--seed",
            "500",
            "--input",
            str(in1),
            str(in2),
            "--output",
            str(out1),
            str(out2),
        ]
    )

    assert out1.read_text() == "@a\nAGC\n+\n123\n@c\nCTG\n+\n***\n"
