import io

from heyfastqlib.argparse_types import (
    GzipFileType,
    set_gzip_output_compresslevel,
)


def test_gzip_file_type_respects_configured_compression_level(monkeypatch, tmp_path):
    captured = {}

    def fake_gzip_open(filename, mode, **kwargs):
        captured["filename"] = filename
        captured["mode"] = mode
        captured["kwargs"] = kwargs
        return io.StringIO()

    monkeypatch.setattr("heyfastqlib.argparse_types.gzip.open", fake_gzip_open)

    output_path = tmp_path / "reads.fastq.gz"
    set_gzip_output_compresslevel(5)
    handle = GzipFileType("w")(str(output_path))
    handle.write("@r1\nACGT\n+\nIIII\n")

    assert captured["filename"] == str(output_path)
    assert captured["mode"] == "wt"
    assert captured["kwargs"]["compresslevel"] == 5


def test_gzip_file_type_uses_default_when_level_not_set(monkeypatch, tmp_path):
    captured = {}

    def fake_gzip_open(filename, mode, **kwargs):
        captured["kwargs"] = kwargs
        return io.StringIO()

    monkeypatch.setattr("heyfastqlib.argparse_types.gzip.open", fake_gzip_open)

    set_gzip_output_compresslevel(None)
    GzipFileType("w")(str(tmp_path / "reads.fastq.gz"))

    assert "compresslevel" not in captured["kwargs"]
