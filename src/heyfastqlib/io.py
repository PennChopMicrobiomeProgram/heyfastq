from .read import R, Read, ReadPair, ReadPipe
from typing import Generator, Iterator, overload, TextIO, Union


def _grouper(iterable: Iterator[str], n: int) -> Iterator[tuple[str, ...]]:
    "Collect data into fixed-length chunks or blocks"
    # grouper('ABCDEFG', 3) --> ABC DEF
    args = [iter(iterable)] * n
    return zip(*args)


def parse_fastq_single(f: TextIO) -> ReadPipe[Read]:
    for desc, seq, _, qual in _grouper(f, 4):
        desc = desc.rstrip()[1:]
        seq = seq.rstrip()
        qual = qual.rstrip()
        yield Read(desc, seq, qual)


def parse_fastq_paired(f1: TextIO, f2: TextIO) -> ReadPipe[ReadPair]:
    yield from zip(parse_fastq_single(f1), parse_fastq_single(f2))


@overload
def parse_fastq(fs: tuple[TextIO]) -> ReadPipe[Read]: ...
@overload
def parse_fastq(fs: tuple[TextIO, TextIO]) -> ReadPipe[ReadPair]: ...


def parse_fastq(
    fs: Union[tuple[TextIO], tuple[TextIO, TextIO]],
) -> Union[ReadPipe[Read], ReadPipe[ReadPair]]:
    # This must stay a plain function (no `yield` here): if it were a
    # generator, `return parse_fastq_single(...)` would end iteration
    # immediately and single-end input would produce no reads.
    if len(fs) == 1:
        return parse_fastq_single(fs[0])
    elif len(fs) == 2:
        return parse_fastq_paired(fs[0], fs[1])
    else:
        raise ValueError("Only single or paired-end FASTQ files are supported.")


def write_fastq_record(f: TextIO, read: Read):
    f.write(f"@{read.desc}\n{read.seq}\n+\n{read.qual}\n")


def write_fastq(
    fs: Union[tuple[TextIO], tuple[TextIO, TextIO]], reads: ReadPipe[R]
) -> None:
    for r in reads:
        if isinstance(r, Read) and len(fs) == 1:
            write_fastq_record(fs[0], r)
        elif isinstance(r, tuple) and len(fs) == 2:
            write_fastq_record(fs[0], r[0])
            write_fastq_record(fs[1], r[1])
        else:
            raise ValueError("Mixing paired/unpaired inputs with files")


def count_reads(f: TextIO) -> int:
    line_count = sum(1 for _ in f)
    return line_count // 4


def parse_seq_ids(f: TextIO) -> Generator[str, None, None]:
    for line in f:
        line = line.strip()
        if line.startswith("#") or (line == ""):
            continue
        seq_id = line.split()[0]
        yield seq_id
