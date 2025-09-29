from types import SimpleNamespace
import io

import pytest

from heyfastqlib.command import (
    ReadStats,
    _iter_with_stats,
    _run_paired_command,
)
from heyfastqlib.paired_reads import map_paired
from heyfastqlib.read import Read


def _make_read(desc: str, seq: str) -> Read:
    return Read(desc, seq, "I" * len(seq))


def _make_fastq_text(reads):
    return "".join(f"@{read.desc}\n{read.seq}\n+\n{read.qual}\n" for read in reads)


def test_read_stats_observes_paired_reads():
    stats = ReadStats()
    paired_read = (
        _make_read("r1/1", "AAAA"),
        _make_read("r1/2", "AAA"),
    )

    stats.observe(paired_read)

    assert stats.total_reads == 2
    assert stats.total_bases == 7
    assert pytest.approx(stats.average_length) == 3.5


def test_iter_with_stats_passthroughs_reads():
    stats = ReadStats()
    pairs = [
        (
            _make_read("r1/1", "AAAA"),
            _make_read("r1/2", "AAA"),
        ),
        (
            _make_read("r2/1", "CC"),
            _make_read("r2/2", "C"),
        ),
    ]

    result = list(_iter_with_stats(iter(pairs), stats))

    assert result == pairs
    assert stats.total_reads == 4
    assert stats.total_bases == 10
    assert pytest.approx(stats.average_length) == 2.5


def test_run_paired_command_logs_input_and_output_stats(capsys):
    mate1_reads = [
        _make_read("r1/1", "AAAA"),
        _make_read("r2/1", "CCCCC"),
    ]
    mate2_reads = [
        _make_read("r1/2", "GGG"),
        _make_read("r2/2", "TTTTTT"),
    ]

    def shorten(read: Read) -> Read:
        return Read(read.desc, read.seq[:-1], read.qual[:-1])

    def transform(reads):
        return map_paired(reads, shorten)

    args = SimpleNamespace(
        input=[
            io.StringIO(_make_fastq_text(mate1_reads)),
            io.StringIO(_make_fastq_text(mate2_reads)),
        ],
        output=[io.StringIO(), io.StringIO()],
    )

    _run_paired_command(args, transform)

    stderr = capsys.readouterr().err.strip().splitlines()
    assert stderr[0] == "Input reads: total=4, average_length=4.50"
    assert stderr[1] == "Output reads: total=4, average_length=3.50"

    output_reads = [
        _make_fastq_text([Read(r.desc, r.seq[:-1], r.qual[:-1]) for r in mate1_reads]),
        _make_fastq_text([Read(r.desc, r.seq[:-1], r.qual[:-1]) for r in mate2_reads]),
    ]

    assert args.output[0].getvalue() == output_reads[0]
    assert args.output[1].getvalue() == output_reads[1]
