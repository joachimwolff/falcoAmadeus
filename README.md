# FalcoAmadeus: a faster, sequencing-aware falco

**FalcoAmadeus** is a performance-focused fork of
[falco](https://github.com/smithlabcode/falco) (itself an emulation of the
popular [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc)
software) that checks large sequencing reads for common problems.

By default it produces **byte-identical output** to upstream falco: a plain
single-end run's `fastqc_data.txt` and `summary.txt` still match the reference
checksums in the upstream test suite (`test/md5sum.txt`). The speed work never
changes results. Everything new below is **opt-in** behind a flag or a new
mode, so existing single-end pipelines are unaffected.

## Why this fork exists

These optimizations were offered to upstream falco as a patch, but the
maintainers declined to merge it. FalcoAmadeus exists to make the faster
implementation available to users who want it, and to add the sequencing-aware
and reporting features below, while keeping the default output identical to the
original.

## What FalcoAmadeus does differently

- **Correct paired-end analysis** (`--paired-with R2`). FastQC and falco only
  process one file at a time, so paired-end PCR-duplicate rates are
  overestimated. FalcoAmadeus analyses the pair together and writes **one
  combined HTML report** with **Combined / Forward / Reverse** tabs, where
  *Combined* pools R1+R2 as a single dataset (GC, adapters, content, quality
  measured over all reads at once, not two independent files). It reports the
  **corrected paired PCR-duplicate rate** — a pair counts as a duplicate only
  when *both* mates are duplicates — and how many R1/R2 reads **actually match**
  into pairs (by read identifier). Per-mate `fastqc_data.txt` files are kept for
  MultiQC, alongside a `paired_duplication.txt` summary.

- **Assay-aware QC** (`--assay <type>`). Tunes each module's pass/warn/fail
  expectations to the library type (`rnaseq`, `amplicon`, `wgbs` bisulfite/EM,
  `chipseq`, `atac`, `smallrna`, `wgs`, …), so genuinely expected base-content
  patterns (e.g. the C→T depletion of bisulfite libraries) are annotated and
  not flagged as failures.

- **A beautified, interactive HTML report.** The report (now named
  `falcoAmadeus_report.html`) has a modern, responsive design: collapsible
  modules with plain-language explanations and FastQC documentation links, a
  **light/dark theme toggle**, **selectable colourblind-safe colour palettes**
  that also recolour the box plot and per-tile heatmap, readable unified
  hover tooltips, and cleaner plots (a proper box-and-whisker per-base-quality
  plot over pass/warn/fail zones, a padded length-distribution histogram, …).

- **MultiQC-friendly.** The `fastqc_data.txt` output is ingested by MultiQC as a
  FastQC report. `--name NAME` sets an explicit sample name (so inputs passed
  with unhelpful names — e.g. Galaxy paired-collection elements all called
  *forward*/*reverse* — stay distinct instead of collapsing), and
  `--paired-data mates|combined|both` chooses whether MultiQC gets the two mates
  separately or a single pooled sample per pair.

- **Multi-threaded gzip parsing** (`-t N`). Gzip decompression runs on a
  background thread overlapped with parsing; output is byte-identical to the
  single-threaded default.

A Galaxy tool wrapper exposing these features is maintained separately as
`tools-iuc/tools/falcoAmadeus`.

## Screenshots

The redesigned report, in light and dark themes:

| Light | Dark |
| :---: | :---: |
| [![Light theme](img/report_light.png)](img/report_light.png) | [![Dark theme](img/report_dark.png)](img/report_dark.png) |

Paired-end combined report (Combined / Forward / Reverse tabs, with the
read-pair summary at the top):

[![Paired-end combined report](img/report_paired.png)](img/report_paired.png)

## Speed improvements

Optimizations are confined to the input parser (`src/StreamReader.*`), which
profiling showed to be the dominant cost:

- Hoist the (virtual) `is_eof()` check out of the per-base loops, where it was
  being called tens of millions of times.
- Cache the per-tile quality/count vectors once per read instead of doing
  three hash-map lookups on every quality base.
- Use a large explicit `setvbuf` buffer for FASTQ/SAM reading, cutting
  `read()` syscalls from ~51,700 to ~220 on a 200&nbsp;MB file; enlarge zlib's
  inflate buffer via `gzbuffer` for gzipped input.
- Single-lookup duplication hashing with buffer reuse instead of a double
  lookup plus a fresh string allocation per read.

Measured best-of-N wall-clock on commodity hardware (page cache warm):

| Input                                   | falco   | FalcoAmadeus | Speedup |
| --------------------------------------- | ------: | -----------: | ------: |
| Uncompressed FASTQ (200&nbsp;MB)        | 0.83&nbsp;s | 0.62&nbsp;s | ~1.34x |
| Diverse FASTQ (118&nbsp;MB, 500k reads) | 1.04&nbsp;s | 0.72&nbsp;s | ~1.44x |
| Gzipped FASTQ (64&nbsp;MB)              | 1.28&nbsp;s | 0.99&nbsp;s | ~1.29x |
| SAM (46&nbsp;MB)                        | 0.22&nbsp;s | 0.18&nbsp;s | ~1.22x |

Peak memory is unchanged: falco already streams input, so resident memory is
constant regardless of file size.

### Parallel (multi-threaded) gzip parsing

With `-t N` (N ≥ 2), gzip decompression runs on a background thread overlapped
with parsing, and per-read statistics are computed in parallel. This helps
**gzipped input only** — uncompressed FASTQ/SAM has no decompression to overlap
— and the output stays **byte-identical to single-threaded** at every thread
count.

Best-of-3 wall-clock on a 675&nbsp;MB gzipped FASTQ (8M reads), 32-core machine:

| Threads (`-t`) | Time | Speedup vs `-t 1` |
| ---: | ---: | ---: |
| 1  | 10.4&nbsp;s | 1.0x |
| 2  |  5.2&nbsp;s | 2.0x |
| 4  |  4.1&nbsp;s | 2.5x |
| 8  |  4.3&nbsp;s | 2.4x |
| 16 |  4.5&nbsp;s | 2.3x |

Speedup peaks around **2.5x at `-t 4`** and then plateaus: a single gzip stream
can only be *inflated* serially, so once decompression is fully overlapped with
parsing there is little left to parallelize, and extra threads only add
scheduling overhead. Combined with the single-threaded parser win above, this
is roughly **3x faster than upstream falco** on gzipped input.

The binary is installed and invoked as `falcoAmadeus`. Apart from the program
name and the added options above, all usage is identical to upstream falco.

Installing FalcoAmadeus
=======================

## Installing through conda

If you use [anaconda](https://anaconda.org) to manage your packages, and the
`conda` binary is in your path, you can install `falcoAmadeus` with

```
$ conda install -c bioconda falcoamadeus
```

The `falcoAmadeus` binary can be found inside the `bin` directory of your
anaconda installer.

> **Note:** the bioconda package is being submitted. Until it lands, build from
> source (below). The bioconda recipe lives in [`recipe/`](recipe/); publishing
> it requires pushing a `v1` release tag to this repository, filling in the
> tarball `sha256` in [`recipe/meta.yaml`](recipe/meta.yaml), and opening a pull
> request against
> [bioconda-recipes](https://github.com/bioconda/bioconda-recipes).

## Installing from a source release

Download a release from the
[releases](https://github.com/joachimwolff/falcoAmadeus/releases) section of
this repository. Extract the `.tar.gz` or `.zip`, `cd` into the extracted
directory (which contains a `configure` script), and run

```
$ ./configure CXXFLAGS="-O3 -Wall"
$ make all
$ make install
```

To install to a specific location, pass `--prefix` to `./configure`:

```
$ ./configure CXXFLAGS="-O3 -Wall" --prefix=/path/to/installation/directory
```

The `falcoAmadeus` binary will be found in the `bin` directory inside the
specified prefix.

## Installing from a cloned repository

For the most recent code (which may contain undocumented bugs), clone the
repository and generate the build system with `autoreconf`:

```
$ git clone https://github.com/joachimwolff/falcoAmadeus.git
$ cd falcoAmadeus
$ ./autogen.sh
$ ./configure CXXFLAGS="-O3 -Wall"
$ make all
$ make install
```

### Required C++ dependencies

[zlib](https://zlib.net) is required to read gzip-compressed FASTQ files. It is
usually installed by default on UNIX systems and is part of the htslib setup,
but can also be installed with apt, brew or conda.

On Ubuntu:

```
$ sudo apt install zlib1g zlib1g-dev
```

### Optional C++ dependencies

[htslib](https://github.com/samtools/htslib) is required to process BAM files.
If not provided, BAM files are treated as unrecognized formats. To build with
htslib, add `--enable-hts`:

```
$ ./configure CXXFLAGS="-O3 -Wall" --enable-hts
```

Running FalcoAmadeus
====================

Run FalcoAmadeus with the following command, where `example.fq` can be replaced
with any FASTQ/FASTQ.gz/SAM/BAM file:

```
$ falcoAmadeus example.fq
```

This generates three files in the output directory:

* `fastqc_data.txt` — a text summary of the QC metrics (ingestible by MultiQC).
* `falcoAmadeus_report.html` — the interactive HTML report.
* `summary.txt` — a tab-separated pass/warn/fail result for each module.

**Paired-end** (one combined report + per-mate data for MultiQC +
`paired_duplication.txt`):

```
$ falcoAmadeus R1.fastq.gz --paired-with R2.fastq.gz -o out/
```

**Assay-aware grading** (here for a bisulfite / EM-seq library):

```
$ falcoAmadeus reads.fastq.gz --assay wgbs -o out/
```

### Options

FalcoAmadeus-specific options:

| Option | Description |
| --- | --- |
| `--paired-with R2` | Analyse the R2 mate together with the R1 positional input; write one combined report and `paired_duplication.txt`. |
| `--paired-data mates\|combined\|both` | With `--paired-with`, which `fastqc_data.txt` files to emit for MultiQC: separate mates (default), a single pooled sample, or both. |
| `--assay TYPE` | Assay-aware pass/warn/fail grading (`rnaseq`, `amplicon`, `wgbs`/`emseq`/`rrbs`, `chipseq`, `atac`, `smallrna`, `wgs`, …; default `generic` is byte-identical to falco). |
| `--conversion-residual F` | For bisulfite assays, the residual unmethylated-base fraction (0–1) above which a possible incomplete-conversion warning is emitted. |
| `--name NAME` | Sample name for the report `Filename` / MultiQC, instead of the input basename. |
| `-t, --threads N` | Overlap gzip decompression with parsing on a background thread (byte-identical output). |

The FastQC-compatible options carried over from falco — `-o/--outdir`,
`-f/--format`, `-c/--contaminants`, `-a/--adapters`, `-l/--limits`, `--nogroup`,
`-s/-subsample`, `-b/-bisulfite`, `-r/-reverse-complement`, the `-skip-*` flags
and `-D`/`-R`/`-S` filename overrides — are unchanged. (The FastQC dummy flags
that falco silently ignored, such as `--casava`, `--nano`, `--extract`,
`--java`, `--kmers`, `--dir` and `--min_length`, have been dropped.)

Run `falcoAmadeus --help` for the complete list.

Copyright and license
=====================

FalcoAmadeus is a **fork** of falco.

* The original falco code is Copyright (C) 2019–2022 Guilherme de Sena Brandine
  and Andrew D. Smith, released under the GNU General Public License. Their
  authorship of the original work is gratefully acknowledged.
* The FalcoAmadeus fork and all modifications to it — the parallel parser, the
  paired-end analysis, assay-aware grading, the redesigned HTML report, the new
  command-line options, and everything else added here — are
  Copyright (C) 2026 Joachim Wolff.

FalcoAmadeus is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

This software is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
