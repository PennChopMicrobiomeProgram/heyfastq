import gzip
import sys
import argparse
from typing import Optional


_gzip_output_compresslevel: Optional[int] = None


def set_gzip_output_compresslevel(level: Optional[int]) -> None:
    """Configure the compression level used for gzip outputs.

    A value of ``None`` falls back to Python's default compression level.
    """

    global _gzip_output_compresslevel
    _gzip_output_compresslevel = level


class GzipFileType(object):
    """Factory for creating optionally gzipped file object types

    Taken from https://github.com/python/cpython/blob/05a370abd6cdfe4b54be60b3b911f3a441026bb2/Lib/argparse.py

    Instances of GzipFileType are typically passed as type= arguments to the
    ArgumentParser add_argument() method.

    Keyword Arguments:
        - mode -- A string indicating how the file is to be opened. Accepts the
            same values as the builtin open() function.
        - bufsize -- The file's desired buffer size. Accepts the same values as
            the builtin open() function.
        - encoding -- The file's encoding. Accepts the same values as the
            builtin open() function.
        - errors -- A string indicating how encoding and decoding errors are to
            be handled. Accepts the same value as the builtin open() function.
    """

    def __init__(self, mode="r", bufsize=-1, encoding=None, errors=None):
        self._mode = mode
        self._bufsize = bufsize
        self._encoding = encoding
        self._errors = errors

    def __call__(self, string):
        # the special argument "-" means sys.std{in,out}
        if string == "-":
            if "r" in self._mode:
                return sys.stdin.buffer if "b" in self._mode else sys.stdin
            elif any(c in self._mode for c in "wax"):
                return sys.stdout.buffer if "b" in self._mode else sys.stdout
            else:
                msg = f'argument "-" with mode {self.mode}'
                raise ValueError(msg)

        # all other arguments are used as file names
        try:
            try:
                with open(string, "rb") as test_f:
                    gzipped = test_f.read(2) == b"\x1f\x8b"
            except FileNotFoundError:
                gzipped = string.endswith(".gz")

            if gzipped:
                gzip_kwargs = {}
                if "w" in self._mode and _gzip_output_compresslevel is not None:
                    gzip_kwargs["compresslevel"] = _gzip_output_compresslevel
                if self._encoding is not None:
                    gzip_kwargs["encoding"] = self._encoding
                if self._errors is not None:
                    gzip_kwargs["errors"] = self._errors

                f = gzip.open(
                    string,
                    f"{self._mode}t",
                    **gzip_kwargs,
                )
            else:
                f = open(
                    string, self._mode, self._bufsize, self._encoding, self._errors
                )

            return f
        except OSError as e:
            args = {"filename": string, "error": e}
            message = f"can't open {args['filename']}: {args['error']}"
            raise argparse.ArgumentTypeError(message % args)

    def __repr__(self):
        args = self._mode, self._bufsize
        kwargs = [("encoding", self._encoding), ("errors", self._errors)]
        args_str = ", ".join(
            [repr(arg) for arg in args if arg != -1]
            + ["%s=%r" % (kw, arg) for kw, arg in kwargs if arg is not None]
        )
        return "%s(%s)" % (type(self).__name__, args_str)


class HFQFormatter(argparse.HelpFormatter):
    def __init__(self, prog, indent_increment=2, max_help_position=25, width=None):
        super().__init__(prog, indent_increment, max_help_position, width)

    # based on ArgumentDefaultsHelpFormatter but with a different search string
    def _get_help_string(self, action):
        help = action.help
        if help is None:
            help = ""

        if "default" not in help:
            if action.default is not argparse.SUPPRESS:
                defaulting_nargs = [argparse.OPTIONAL, argparse.ZERO_OR_MORE]
                if action.option_strings or action.nargs in defaulting_nargs:
                    help += " (default: %(default)s)"
        return help
