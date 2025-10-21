import gzip
from contextlib import contextmanager
from pathlib import Path

from heyfastqlib.argparse_types import GzipFileType


@contextmanager
def gzft_open(gzft: GzipFileType, path: Path):
    handle, closer = gzft(str(path))
    try:
        yield handle
    finally:
        if closer:
            closer()
        handle.close()


def test_gzipfiletype_init(tmp_path):
    gzft = GzipFileType()
    assert gzft._mode == "r"


def test_gzipfiletype_call(tmp_path):
    gzft = GzipFileType()
    source = tmp_path / "test_in.txt"
    source.write_text("test")
    with gzft_open(gzft, source) as handle:
        assert handle.read() == "test"

    gzftw = GzipFileType(mode="w")
    dest = tmp_path / "test_out.txt"
    with gzft_open(gzftw, dest) as handle:
        handle.write("test")
    with gzft_open(gzft, dest) as handle:
        assert handle.read() == "test"


def test_gzipfiletype_call_gz(tmp_path):
    gzft = GzipFileType()
    gz_source = tmp_path / "test_in.txt.gz"
    with gzip.open(gz_source, "wt") as handle:
        handle.write("test")
    with gzft_open(gzft, gz_source) as handle:
        assert handle.read() == "test"

    gzftw = GzipFileType(mode="wt")
    gz_dest = tmp_path / "test_out.txt.gz"
    with gzft_open(gzftw, gz_dest) as handle:
        handle.write("test")
    with gzft_open(gzft, gz_dest) as handle:
        assert handle.read() == "test"
