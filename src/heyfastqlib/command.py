import argparse
import operator
import signal
import sys

from . import __version__
from .util import (
    subsample,
)
from .io import (
    parse_fastq_paired,
    write_fastq_paired,
    parse_seq_ids,
)
from .paired_reads import map_paired, filter_paired
from .read import (
    trim,
    kscore_ok,
    length_ok,
    seq_id_ok,
    trim_moving_average,
    trim_ends,
)
from .argparse_types import GzipFileType, HFQFormatter


def subsample_subcommand(args):
    reads = parse_fastq_paired(args.input)
    out_reads = subsample(reads, args.n)
    write_fastq_paired(args.output, out_reads)


def trim_fixed_subcommand(args):
    reads = parse_fastq_paired(args.input)
    out_reads = map_paired(reads, trim, end_idx=args.length)
    write_fastq_paired(args.output, out_reads)


def trim_qual_subcommand(args):
    reads = parse_fastq_paired(args.input)
    trimmed_moving_average_reads = map_paired(
        reads, trim_moving_average, k=args.window_width, threshold=args.window_threshold
    )
    trimmed_ends_reads = map_paired(
        trimmed_moving_average_reads,
        trim_ends,
        threshold_start=args.start_threshold,
        threshold_end=args.end_threshold,
    )
    filtered_reads = filter_paired(
        trimmed_ends_reads, length_ok, threshold=args.min_length
    )
    write_fastq_paired(args.output, filtered_reads)


def filter_length_subcommand(args):
    reads = parse_fastq_paired(args.input)
    cmp = operator.lt if args.less else operator.ge
    out_reads = filter_paired(reads, length_ok, threshold=args.length, cmp=cmp)
    write_fastq_paired(args.output, out_reads)


def filter_kscore_subcommand(args):
    reads = parse_fastq_paired(args.input)
    out_reads = filter_paired(
        reads, kscore_ok, k=args.kmer_size, min_kscore=args.min_kscore
    )
    write_fastq_paired(args.output, out_reads)


def filter_seq_ids_subcommand(args):
    seq_ids = set(parse_seq_ids(args.idsfile))
    reads = parse_fastq_paired(args.input)
    out_reads = filter_paired(reads, seq_id_ok, seq_ids=seq_ids, keep=args.keep_ids)
    write_fastq_paired(args.output, out_reads)


fastq_io_parser = argparse.ArgumentParser(add_help=False, formatter_class=HFQFormatter)
fastq_io_parser.add_argument(
    "--input",
    type=GzipFileType("r"),
    nargs="*",
    default=[sys.stdin],
    help="Input FASTQs, can be gzipped (default: stdin)",
)
fastq_io_parser.add_argument(
    "--output",
    type=GzipFileType("w"),
    nargs="*",
    default=[sys.stdout],
    help="Output FASTQs, can be gzipped (default: stdout)",
)


def heyfastq_main(argv=None):
    # Ignore SIG_PIPE and don't throw exceptions on it
    # newbebweb.blogspot.com/2012/02/python-head-ioerror-errno-32-broken.html
    # Try/catch to not fail on Windows
    # https://github.com/t2mune/mrtparse/issues/18
    try:
        signal.signal(signal.SIGPIPE, signal.SIG_DFL)
    except AttributeError:
        pass

    main_parser = argparse.ArgumentParser()
    main_parser.add_argument(
        "-v",
        "--version",
        action="version",
        version=f"{__version__}",
    )
    subparsers = main_parser.add_subparsers(title="Subcommands", required=True)

    trim_fixed_parser = subparsers.add_parser(
        "trim-fixed",
        parents=[fastq_io_parser],
        formatter_class=HFQFormatter,
        help="Trim reads to fixed length",
    )
    trim_fixed_parser.add_argument(
        "--length",
        type=int,
        default=100,
        help="Length of output reads",
    )
    trim_fixed_parser.set_defaults(func=trim_fixed_subcommand)

    trim_qual_parser = subparsers.add_parser(
        "trim-qual",
        parents=[fastq_io_parser],
        formatter_class=HFQFormatter,
        help="Trim reads based on quality scores",
    )
    trim_qual_parser.add_argument(
        "--window-width", type=int, default=4, help="Sliding window width"
    )
    trim_qual_parser.add_argument(
        "--window-threshold",
        type=float,
        default=15,
        help="Sliding window mean quality threshold",
    )
    trim_qual_parser.add_argument(
        "--start-threshold",
        type=float,
        default=3,
        help="Quality threshold for trimming start of read",
    )
    trim_qual_parser.add_argument(
        "--end-threshold",
        type=float,
        default=3,
        help="Quality threshold for trimming end of read",
    )
    trim_qual_parser.add_argument(
        "--min-length",
        type=int,
        default=36,
        help="Minimum length after quality trimming",
    )
    trim_qual_parser.set_defaults(func=trim_qual_subcommand)

    filter_length_parser = subparsers.add_parser(
        "filter-length",
        parents=[fastq_io_parser],
        formatter_class=HFQFormatter,
        help="Filter reads by length",
    )
    filter_length_parser.add_argument(
        "--length",
        type=int,
        default=100,
        help="Length threshold",
    )
    filter_length_parser.add_argument(
        "--less",
        action="store_true",
        help=(
            "Keep reads that are less than the specified length "
            "(default: keep greater than or equal to length)"
        ),
    )
    filter_length_parser.set_defaults(func=filter_length_subcommand)

    filter_kscore_parser = subparsers.add_parser(
        "filter-kscore",
        parents=[fastq_io_parser],
        formatter_class=HFQFormatter,
        help="Filter reads by komplexity score",
    )
    filter_kscore_parser.add_argument(
        "--kmer-size", type=int, default=4, help="Kmer size"
    )
    filter_kscore_parser.add_argument(
        "--min-kscore",
        type=float,
        default=0.55,
        help="Minimum komplexity score",
    )
    filter_kscore_parser.set_defaults(func=filter_kscore_subcommand)

    filter_seq_ids_parser = subparsers.add_parser(
        "filter-seqids",
        parents=[fastq_io_parser],
        formatter_class=HFQFormatter,
        help="Filter reads by sequence id",
    )
    filter_seq_ids_parser.add_argument(
        "idsfile",
        type=argparse.FileType("r"),
        help="File containing sequence ids, one per line",
    )
    filter_seq_ids_parser.add_argument(
        "--keep-ids",
        action="store_true",
        help="Keep, rather than remove ids in list",
    )
    filter_seq_ids_parser.set_defaults(func=filter_seq_ids_subcommand)

    subsample_parser = subparsers.add_parser(
        "subsample",
        parents=[fastq_io_parser],
        formatter_class=HFQFormatter,
        help="Select random reads",
    )
    subsample_parser.add_argument("--n", type=int, default=1000, help="Number of reads")
    subsample_parser.set_defaults(func=subsample_subcommand)

    args = main_parser.parse_args(argv)
    if args.input is None:  # pragma: no cover
        args.input = sys.stdin
    if args.output is None:  # pragma: no cover
        args.output = sys.stdout
    args.func(args)
