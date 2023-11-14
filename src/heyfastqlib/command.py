import argparse
import signal
import sys

from .util import (
    subsample,
)
from .io import (
    parse_fastq_paired, write_fastq_paired,
)
from .paired_reads import (
    map_paired, filter_paired, trim_fixed, kscore_ok,
    )

def subsample_subcommand(args):
    reads = parse_fastq_paired(args.input)
    out_reads = subsample(reads, args.n)
    write_fastq_paired(args.output, out_reads)

def trim_fixed_subcommand(args):
    reads = parse_fastq_paired(args.input)
    out_reads = map_paired(reads, trim_fixed, length=args.length)
    write_fastq_paired(args.output, out_reads)

def filter_kscore_subcommand(args):
    reads = parse_fastq_paired(args.input)
    out_reads = filter_paired(
        reads, kscore_ok, k=args.kmer_size, min_kscore=args.min_kscore)
    write_fastq_paired(args.output, out_reads)

fastq_io_parser = argparse.ArgumentParser(add_help=False)
fastq_io_parser.add_argument(
    "--input", type=argparse.FileType('r'), nargs="*", default=[sys.stdin],
    help=(
        "Input FASTQ, can specify more than one file for paired reads "
        "(default: stdin)"))
fastq_io_parser.add_argument(
    "--output", type=argparse.FileType('w'), nargs="*", default=[sys.stdout],
    help=(
        "Output FASTQ, can specify more than one file for paired reads "
        "(default: stdout)"))

def heyfastq_main(argv=None):
    # Ignore SIG_PIPE and don't throw exceptions on it
    # newbebweb.blogspot.com/2012/02/python-head-ioerror-errno-32-broken.html
    signal.signal(signal.SIGPIPE, signal.SIG_DFL)

    main_parser = argparse.ArgumentParser()
    subparsers = main_parser.add_subparsers(
        title="Subcommands", required=True)

    trim_fixed_parser = subparsers.add_parser(
        "trim-fixed", parents=[fastq_io_parser],
        help="Trim sequences to fixed length")
    trim_fixed_parser.add_argument(
        "--length", type=int, default=100,
        help="Length of output sequences")
    trim_fixed_parser.set_defaults(func=trim_fixed_subcommand)

    filter_kscore_parser = subparsers.add_parser(
        "filter-kscore", parents=[fastq_io_parser],
        help="Filter read pairs by komplexity score")
    filter_kscore_parser.add_argument(
        "--kmer-size", type=int, default=4,
        help="Kmer size (default: %(default)s)")
    filter_kscore_parser.add_argument(
        "--min-kscore", type=float, default=0.55,
        help="Minimum komplexity score (default: %(default)s)")
    filter_kscore_parser.set_defaults(func=filter_kscore_subcommand)
    
    subsample_parser = subparsers.add_parser(
        "subsample", parents=[fastq_io_parser],
        help='Select random sequences')
    subsample_parser.add_argument(
        "--n", type=int, default=1000,
        help="Number of sequences (default: %(default)s)")
    subsample_parser.set_defaults(func=subsample_subcommand)

    args = main_parser.parse_args(argv)
    if args.input is None: # pragma: no cover
        args.input = sys.stdin
    if args.output is None: # pragma: no cover
        args.output = sys.stdout
    args.func(args)
