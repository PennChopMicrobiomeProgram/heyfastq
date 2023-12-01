import gzip
from heyfastqlib.argparse_types import GzipFileType

def test_gzipfiletype_init(tmp_path):
    gzft = GzipFileType()
    assert gzft._mode == 'r'

def test_gzipfiletype_call(tmp_path):
    gzft = GzipFileType()
    with open(tmp_path / "test_in.txt", "w") as f:
        f.write("test")
    with gzft(str(tmp_path / "test_in.txt")) as f:
        assert f.read() == "test"
    
    gzftw = GzipFileType(mode = "w")
    with gzftw(str(tmp_path / "test_out.txt")) as f:
        f.write("test")
    with gzft(str(tmp_path / "test_out.txt")) as f:
        assert f.read() == "test"

def test_gzipfiletype_call_gz(tmp_path):
    gzft = GzipFileType()
    with gzip.open(tmp_path / "test_in.txt.gz", "wt") as f:
        f.write("test")
    with gzft(str(tmp_path / "test_in.txt.gz")) as f:
        assert f.read() == "test"
    
    gzftw = GzipFileType(mode = "wt")
    with gzftw(str(tmp_path / "test_out.txt.gz")) as f:
        f.write("test")
    with gzft(str(tmp_path / "test_out.txt.gz")) as f:
        assert f.read() == "test"