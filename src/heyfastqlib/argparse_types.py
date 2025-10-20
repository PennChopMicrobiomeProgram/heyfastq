import argparse
import gzip
import io
import os
import shutil
import subprocess
import sys

from typing import Iterable, Optional


_DEFAULT_THREAD_COUNT = str(min(32, os.cpu_count() or 32))
# TEST
_DEFAULT_THREAD_COUNT = os.environ.get("HFQ_BENCH_PIGZ_THREADS", None)


def _normalise_gzip_mode(mode: str) -> str:
    if "t" in mode or "b" in mode:
        return mode
    return f"{mode}t"


def _gzip_reader_process(path: str) -> Optional[subprocess.Popen]:
    cmd = None

    pigz = shutil.which("pigz")
    if pigz:
        cmd = [pigz, "-p", _DEFAULT_THREAD_COUNT, "-cd", path]
    else:
        gzip_bin = shutil.which("gzip")
        if gzip_bin:
            cmd = [gzip_bin, "-cd", path]


    print(cmd)
    if cmd is None:
        return None

    return subprocess.Popen(cmd, stdout=subprocess.PIPE)


def _gzip_writer_process(destination: str) -> Optional[tuple[subprocess.Popen, Iterable[io.IOBase]]]:
    cmd = None

    pigz = shutil.which("pigz")
    if pigz:
        cmd = [pigz, "-p", _DEFAULT_THREAD_COUNT, "-c"]
    else:
        gzip_bin = shutil.which("gzip")
        if gzip_bin:
            cmd = [gzip_bin, "-c"]

    print(cmd)
    if cmd is None:
        return None

    output_handle = open(destination, "wb")
    proc = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=output_handle)
    return proc, (output_handle,)


class _SubprocessFileWrapper:
    def __init__(self, stream, process: Optional[subprocess.Popen], extra_closers: Iterable[io.IOBase] = ()): 
        self._stream = stream
        self._process = process
        self._extra_closers = list(extra_closers)
        self._closed = False

    def __getattr__(self, item):
        return getattr(self._stream, item)

    def __iter__(self):
        return iter(self._stream)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, tb):
        self.close()
        return False

    @property
    def closed(self):
        if self._closed:
            return True
        return getattr(self._stream, "closed", False)

    def close(self):
        if self._closed:
            return

        stream_exc = None
        try:
            self._stream.close()
        except Exception as exc:  # pragma: no cover - propagate original exception
            stream_exc = exc

        closer_exc = None
        for closer in self._extra_closers:
            try:
                closer.close()
            except Exception as exc:  # pragma: no cover - propagate original exception
                if closer_exc is None:
                    closer_exc = exc

        process_exc = None
        if self._process is not None:
            return_code = self._process.wait()
            if return_code != 0:
                process_exc = subprocess.CalledProcessError(return_code, self._process.args)

        self._closed = True

        if stream_exc is not None:
            raise stream_exc
        if closer_exc is not None:
            raise closer_exc
        if process_exc is not None:
            raise process_exc


def _wrap_reader(process: Optional[subprocess.Popen], mode: str, encoding, errors):
    if process is None or process.stdout is None:
        return None

    stream = process.stdout
    if "b" not in mode:
        stream = io.TextIOWrapper(stream, encoding=encoding, errors=errors)

    return _SubprocessFileWrapper(stream, process)


def _wrap_writer(process_info: Optional[tuple[subprocess.Popen, Iterable[io.IOBase]]], mode: str, encoding, errors):
    if process_info is None:
        return None

    process, extra_closers = process_info
    if process.stdin is None:
        for closer in extra_closers:
            closer.close()
        process.terminate()
        raise RuntimeError("Failed to initialise compression process")

    stream = process.stdin
    if "b" not in mode:
        stream = io.TextIOWrapper(stream, encoding=encoding, errors=errors)

    return _SubprocessFileWrapper(stream, process, extra_closers)


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
                if "r" in self._mode:
                    wrapper = _wrap_reader(
                        _gzip_reader_process(string),
                        self._mode,
                        self._encoding,
                        self._errors,
                    )
                    if wrapper is None:
                        f = gzip.open(
                            string,
                            _normalise_gzip_mode(self._mode),
                            self._bufsize,
                            self._encoding,
                            self._errors,
                        )
                    else:
                        f = wrapper
                elif any(c in self._mode for c in "wax"):
                    wrapper = _wrap_writer(
                        _gzip_writer_process(string),
                        self._mode,
                        self._encoding,
                        self._errors,
                    )
                    if wrapper is None:
                        f = gzip.open(
                            string,
                            _normalise_gzip_mode(self._mode),
                            self._bufsize,
                            self._encoding,
                            self._errors,
                        )
                    else:
                        f = wrapper
                else:
                    msg = f'argument "{string}" with mode {self._mode}'
                    raise ValueError(msg)
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
