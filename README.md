# heyfastq

<!-- Badges start -->
[![CI](https://github.com/PennChopMicrobiomeProgram/heyfastq/actions/workflows/ci.yml/badge.svg)](https://github.com/PennChopMicrobiomeProgram/heyfastq/actions/workflows/ci.yml)
[![PyPI version](https://badge.fury.io/py/heyfastq.svg)](https://pypi.org/project/heyfastq/)
[![DockerHub](https://img.shields.io/docker/pulls/chopmicrobiome/heyfastq)](https://hub.docker.com/repository/docker/chopmicrobiome/heyfastq/)
<!-- Badges end -->

FASTQ sequence file utilities, written in pure Python, with no
dependencies.

## Summary

The package comes with one program, `heyfastq`, which provides
utilities for single or paired FASTQ files.

## Installation

Install from PyPi with:

```bash
pip install heyfastq
```

Or get the dev version from GitHub:

```bash
git clone https://github.com/kylebittinger/heyfastq.git
pip install .
```

## Usage

Run `heyfastq -h` to learn more about usage options.

### As dependency

To include `heyfastqlib` in your own Python code, include it as a dependency (only available via GH atm). Then import the objects/functions you need and enjoy your fastqs.

e.g. your very own wrapper for fastq_parse that only returns the sequences

```
import gzip
from heyfastqlib.io import parse_fastq
from typing import Generator

def fastq_reads(path: str) -> Generator[str, None, None]:
    with gzip.open("/path/to/my_file_R1.fastq.gz", "rt") as f:
        return (read.seq for read in parse_fastq(f.readlines()))
```
