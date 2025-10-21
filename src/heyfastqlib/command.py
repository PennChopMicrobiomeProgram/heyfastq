import argparse
import operator
import signal
import sys
from . import __version__
from .argparse_types import GzipFileType, HFQFormatter
from .io import (
    count_reads,
    parse_fastq,
    write_fastq,
    parse_seq_ids,
)
from .pipelines import filter_reads, map_reads, subsample_reads
from .read import (
    trim,
    kscore_ok,
    length_ok,
    seq_id_ok,
    trim_moving_average,
    trim_ends,
)


def subsample_subcommand(args):
    num_reads = count_reads(args.input[0])
    args.input[0].seek(0)
    write_fastq(
        args.output,
        subsample_reads(parse_fastq(args.input), num_reads, args.n, args.seed),
    )


def trim_fixed_subcommand(args):
    counter = {"input_reads": 0, "input_bases": 0, "output_reads": 0, "output_bases": 0}
    write_fastq(
        args.output,
        map_reads(
            parse_fastq(args.input),
            trim,
            counter,
            start_idx=0,
            end_idx=args.length,
            threads=args.threads,
            chunk_size=args.chunk_size,
        ),
    )


def trim_qual_subcommand(args):
    length_counter = {
        "input_reads": 0,
        "input_bases": 0,
        "output_reads": 0,
        "output_bases": 0,
    }
    trim_ends_counter = {
        "input_reads": 0,
        "input_bases": 0,
        "output_reads": 0,
        "output_bases": 0,
    }
    trim_avg_counter = {
        "input_reads": 0,
        "input_bases": 0,
        "output_reads": 0,
        "output_bases": 0,
    }
    reads = parse_fastq(args.input)
    reads = map_reads(
        reads,
        trim_moving_average,
        trim_avg_counter,
        k=args.window_width,
        threshold=args.window_threshold,
        threads=args.threads,
        chunk_size=args.chunk_size,
    )
    reads = map_reads(
        reads,
        trim_ends,
        trim_ends_counter,
        threshold_start=args.start_threshold,
        threshold_end=args.end_threshold,
        threads=args.threads,
        chunk_size=args.chunk_size,
    )
    reads = filter_reads(
        reads,
        length_ok,
        length_counter,
        threshold=args.min_length,
        threads=args.threads,
        chunk_size=args.chunk_size,
    )
    write_fastq(args.output, reads)


def filter_length_subcommand(args):
    cmp = operator.lt if args.less else operator.ge
    counter = {"input_reads": 0, "input_bases": 0, "output_reads": 0, "output_bases": 0}
    write_fastq(
        args.output,
        filter_reads(
            parse_fastq(args.input),
            length_ok,
            counter,
            threshold=args.length,
            cmp=cmp,
            threads=args.threads,
            chunk_size=args.chunk_size,
        ),
    )


def filter_kscore_subcommand(args):
    counter = {"input_reads": 0, "input_bases": 0, "output_reads": 0, "output_bases": 0}
    write_fastq(
        args.output,
        filter_reads(
            parse_fastq(args.input),
            kscore_ok,
            counter,
            k=args.kmer_size,
            min_kscore=args.min_kscore,
            threads=args.threads,
            chunk_size=args.chunk_size,
        ),
    )


def filter_seq_ids_subcommand(args):
    seq_ids = set(parse_seq_ids(args.idsfile))
    counter = {"input_reads": 0, "input_bases": 0, "output_reads": 0, "output_bases": 0}
    write_fastq(
        args.output,
        filter_reads(
            parse_fastq(args.input),
            seq_id_ok,
            counter,
            seq_ids=seq_ids,
            keep=args.keep_ids,
            threads=args.threads,
            chunk_size=args.chunk_size,
        ),
    )


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
fastq_io_parser.add_argument(
    "--threads", type=int, default=1, help="Number of threads to use (default: 1)"
)
fastq_io_parser.add_argument(
    "--chunk-size",
    type=int,
    default=1000,
    help="Number of reads processed per worker chunk (default: 1000)",
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
    subsample_parser.add_argument("--seed", type=int, help="Random seed")
    subsample_parser.set_defaults(func=subsample_subcommand)

    args = main_parser.parse_args(argv)

    # This closers list is a pretty convoluted mechanism to ensure that all opened files and pipes are closed after use
    # It handles everything from sys.stdin/out (no closing necessary) to subprocess pipes (need to close stream handlers and wait for process to end)
    # So we attach a closer function to each opened input/output file handler and call them all at the end
    closers = []
    if args.input is None:
        args.input = sys.stdin
    else:
        for i in args.input:
            closers.append(i[1])
        args.input = [i[0] for i in args.input]
    if args.output is None:
        args.output = sys.stdout
    else:
        for o in args.output:
            closers.append(o[1])
        args.output = [o[0] for o in args.output]
    if args.threads is None:
        args.threads = 1
    args.func(args)

    for c in closers:
        c()
