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

## Dev

Heyfastq is built around the idea of piping reads (or read pairs) through filter and map functions. The fundamental unit that moves through heyfastq pipelines is the `R` object, which can be either a `Read` or a `ReadPair`. These generic `R`s move through functions that take in `ReadPipe`s and output `ReadPipe`s, allowing for easy composition of pipelines. 

```
from heyfastqlib.read import Read, ReadPair, R, ReadPipe
from heyfastqlib.pipelines import filter_reads, map_reads

def unit_filter(r: R) -> bool:
  return True

def unit_map(r: R) -> R:
  return r

input_fastq = (r for r in [Read("1", "ACTG", "HHHH"), Read("2", "GTCA", "HHHH"), Read("3", "AAAA", "####")])
filter_counter = {"input_reads": 0, "input_bases": 0, "output_reads": 0, "output_bases": 0}
map_counter = {"input_reads": 0, "input_bases": 0, "output_reads": 0, "output_bases": 0}

output_fastq = map_reads(filter_reads(input_fastq, unit_filter, filter_counter), unit_map, map_counter)
```

This is all well and good, but how do we actually deal with fastq files, not just objects already in python?

```
from heyfastqlib.io import parse_fastq, write_fastq

with open("r1.fq") as f_in, open("o1.fq", "w") as f_out:
  write_fastq(f_out, parse_fastq(f_in))

with open("r1.fq") as f1_in, open("r2.fq") as f2_in, open("o1.fq", "w") as f1_out, open("o2.fq", "w") as f2_out:
  write_fastq((f1_out, f2_out), parse_fastq((f1_in, f2_in)))
```

Putting it all together, assuming we've imported/created the objects from above, let's read our fastqs then apply our filter then our map and finally write them to fastqs again:

```
with open("r1.fq") as f1_in, open("r2.fq") as f2_in, open("o1.fq", "w") as f1_out, open("o2.fq", "w") as f2_out:
  write_fastq((f1_out, f2_out), map_reads(filter_reads(parse_fastq((f1_in, f2_in)), unit_filter, filter_counter), unit_map, map_counter))
```