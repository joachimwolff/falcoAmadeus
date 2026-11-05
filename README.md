# FalcoAmadeus: a faster, sequencing-aware falco

**FalcoAmadeus** is a performance-focused fork of
[falco](https://github.com/smithlabcode/falco) (itself an emulation of the
popular [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc)
software) that checks large sequencing reads for common problems.

By default it produces the same output as upstream falco 1.3.2 with one
corrected grade: Per base N content is graded as in FastQC (warn above 5 %,
fail above 20 % N at a position), where falco 1.x never warns or fails because
it compares a fraction with the percent limits. With
`--legacy-n-content-grade`, a plain single-end run's `fastqc_data.txt` and
`summary.txt` are byte-identical to those of falco 1.3.2 (checked in CI). The
speed work never changes results. Everything else new below is **opt-in**
behind a flag or a new mode.

## Why this fork exists

These optimizations were offered to upstream falco as a patch, but the
maintainers declined to merge it. FalcoAmadeus exists to make the faster
implementation available to users who want it, and to add the sequencing-aware
and reporting features below, while keeping the output identical to the
original apart from the corrected N-content grade.

## What FalcoAmadeus does differently

- **Paired-end analysis** (`--paired-with R2`). FastQC and falco process one
  file at a time, so sequence duplication is assessed per mate: two fragments
  that share only one end count as duplicates. FalcoAmadeus analyses the pair
  together and writes **one combined HTML report** with **Combined / Forward /
  Reverse** tabs, where *Combined* pools R1+R2 as a single dataset (GC,
  adapters, content, quality measured over all reads at once). It first checks
  that the mates are in matching order (identical read identifiers) and
  rejects inputs where they are not (`--allow-unmatched-pairs` overrides).
  `paired_duplication.txt` reports the **paired sequence duplicate rate**,
  where a pair is a copy only when *both* mates match, counted exactly over
  all reads and therefore independent of read order, next to the per-mate
  rates and the FastQC-style estimates for comparison with FastQC/falco.
  On request the report adds an **error- and UMI-aware rate**:
  `--dup-mismatches N` lets the mates of two copies differ in up to N bases in
  total (copies with sequencing errors), and `--umi-r1`/`--umi-r2 START:LEN`
  or `--umi-in-name` keep pairs with different unique molecular identifiers
  apart (`--umi-mismatches` tolerates errors in the UMI). Candidates are found
  through k-mer blocks shared by both mates and confirmed base by base, so the
  rate stays exact for the stated tolerance and independent of read order.
  The same options work for single-end input, where a read is the unit and
  the rate is written to `aware_duplication.txt`.

- **Optical duplicates.** `--optical` counts the duplicates that are
  copies of one cluster read twice by the camera, with Picard's rule: copies
  on the same tile at most `--optical-distance` pixels apart (100 as Picard;
  2500 for patterned flow cells), one per cluster of close copies kept, the
  location from the tile:x:y of the read 1 name, the two mate orders apart,
  x and y kept in 16 bits as Picard's MarkDuplicates does. The sets are the
  groups of copies of the aware pass (exact, or with `--dup-mismatches` and
  the UMI options), so no alignment is needed; pairs that would not align
  (an insert under 30 bases, a mate mostly N) are left out, as Picard sees
  aligned pairs only. Where these sets equal Picard's, the count equals
  Picard's; on public runs see `paper/data/optical_validation`.

- **UMI quality.** With a UMI option (`--umi-r1`, `--umi-r2`,
  `--umi-in-name`) the run also writes `umi_qc.txt`: UMIs with N,
  homopolymer UMIs (poly-G from dark cycles), distinct UMIs and the share
  of the possible sequences in use, the most frequent UMI and the base
  composition per UMI position; the families (groups of copies with the
  same UMI and sequence within `--dup-mismatches` and `--umi-mismatches`)
  with their sizes, the reads in families of 2 and 3 or more, and the
  families per distinct UMI (high values: molecules share UMIs); and the
  base quality of in-read UMIs (first million reads).

- **Machine-readable summary.** Every run writes
  `falcoAmadeus_summary.json` (all metrics per view and of the falcoAmadeus
  reports, plus the share of bases with Phred quality 20 and 30 or more) and
  `falcoAmadeus_mqc.json`, which MultiQC shows as general statistics columns
  (Q30 bases, paired and aware duplication, complexity, insert size, adapter
  read-through and dimers, screen, optical duplicates, UMI family size and
  UMI space used, Q20 trimming).

- **Acceptance rules.** `--acceptance rules.tsv` grades a sample against
  lab limits (lines `metric >= warn fail` or `metric <= warn fail`, the
  metric a path into the JSON summary such as
  `screen.references.human.found_pct`, or a MultiQC column such as
  `q30_bases_pct`) and writes `acceptance.txt` with each rule and an overall
  PASS, WARN or FAIL; a metric the run did not produce counts as FAIL.
  `--acceptance-exit` makes a failed sample end the run with status 2.

- **Contamination screen.** `--screen DB` estimates which share of the
  reads comes from each genome of a reference panel (human, mouse, rat,
  fly, worm, yeast, Arabidopsis, E. coli, PhiX, lambda, mitochondria, rRNA,
  mycoplasma, vectors) without alignment: canonical 31-mers selected by
  hash (FracMinHash), low-complexity k-mers skipped on both sides, in
  `screen.txt` (found per reference, found only there, no reference).
  `--screen-build DB --screen-refs name=fasta,...` builds the database from
  any FASTA files (about 40 MB and 15 s for the panel above);
  `paper/scripts/screen_validation` lists the public references used.

- **Trimming preview.** `--trim-preview` writes `trimming_preview.txt`:
  what cutadapt's 3' quality trimming (`-q` at 10, 20 and 30) would remove
  and how many reads would end up shorter than 20, 36 and 50 bases (the
  same lengths as cutadapt; one section per mate for paired input), and an
  estimate of adapter trimming with the adapter of the list found most
  (`paper/data/trim_preview_validation`).

- **Insert size and adapter read-through (paired).** From the overlap of
  read 1 with the reverse complement of read 2 in the first million pairs,
  `paired_insert_size.txt` gives the insert-size histogram (as fastp's
  overlap analysis; equal histograms on the same pairs), the share of pairs
  whose insert is shorter than a read, the adapter each read continues into
  (named when it is a known adapter), and the adapter dimers.
  `--insert-size-pairs N` sets the number of pairs (0: all).

- **Library complexity.** From how many distinct sequences were seen once,
  twice and so on, `paired_duplication.txt` reports how many distinct pairs
  to expect at 2x to 100x the sequenced pairs, so a library can be judged
  for further sequencing before alignment. The estimate is a port of the
  `preseq lc_extrap` yield curve (quick mode) and gives the same values as
  preseq on the same counts. `--complexity` adds it for single-end input
  (`complexity.txt`).
  Sequence duplication reflects PCR duplicates as well as natural redundancy
  (deep coverage, highly expressed or enriched loci). Per-mate
  `fastqc_data.txt` files are kept for MultiQC.

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

- **Multi-threaded parsing of every input format** (`-t N`). FASTQ, gzipped
  FASTQ, SAM and BAM all parse in parallel, and the same fast drivers serve
  single-threaded runs too. Output is byte-identical to upstream falco at
  every thread count.

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

FalcoAmadeus replaces falco's serial input path with parallel, format-aware
drivers, while keeping the output byte-identical (verified against the
upstream checksum suite and against single-threaded output at every thread
count, for every change):

- **Parallel parsing for every input format.** Uncompressed FASTQ and SAM are
  memory-mapped and parsed in place across record-aligned byte ranges. BAM is
  decompressed and parsed in overlapped, double-buffered chunks (BGZF blocks
  inflate in parallel with libdeflate). Plain gzip inflates with ISA-L into
  bounded chunks that are parsed in parallel while the next chunk
  decompresses. Single-core runs use the same drivers, so `-t 1` is fast too.
- **SIMD statistics kernels** (AVX-512 on x86, NEON on ARM, scalar fallback):
  the sequence pass (separator scan, N detection, GC, per-position base
  composition into small cache-resident accumulator tiles), the quality pass,
  chain-free adapter window probes, and BAM's packed-base decode including
  reverse-complement.
- **Order-exact parallel bookkeeping.** Everything order-dependent (the
  duplication tracker's first-100k-unique semantics, the floating-point GC
  model, read subsampling) reproduces the serial trajectory exactly; that is
  what keeps the output byte-identical rather than approximately equal.
- **Bounded memory.** Inputs stream through reused chunk buffers with
  drop-behind page management, so peak memory stays flat (roughly 20 to
  200&nbsp;MB for single-end input and under 500&nbsp;MB in paired mode,
  depending on format and thread count) no matter how large the input is.

### Benchmark against falco 2.0.3, falco 1.3.2 and FastQC 0.12.1

AMD Ryzen Threadripper 7970X (32 cores, 64 threads), warm page cache, whole
machine; medians of two (single-end) or three (paired-end) runs. falco 2.0.3
was run with its default duplication and with `--orig-dups` (the FastQC
estimate, as falcoAmadeus reports it); the faster of the two is shown. falco
1.3.2 and FastQC run single-threaded. Wall time, with peak memory in
parentheses:

**Single-end, 29.9 million reads of 1 to 76 bp (SRR14777166)**

| Input | FastQC `-t 1` | falco 1.3.2 | falco 2.0.3 `-t 1` | falcoAmadeus `-t 1` | falco 2.0.3 `-t 8` | falcoAmadeus `-t 8` | falco 2.0.3 `-t 32` | falcoAmadeus `-t 32` |
| ----- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| FASTQ      | 32.6&nbsp;s | 20.4&nbsp;s | 5.70&nbsp;s (287&nbsp;MB) | **3.59&nbsp;s** (31&nbsp;MB) | 1.04&nbsp;s (299&nbsp;MB) | **0.57&nbsp;s** (95&nbsp;MB) | 0.49&nbsp;s (322&nbsp;MB) | **0.23&nbsp;s** (103&nbsp;MB) |
| gzip FASTQ | 34.6&nbsp;s | 23.2&nbsp;s | 9.65&nbsp;s (298&nbsp;MB) | **4.35&nbsp;s** (63&nbsp;MB) | 5.28&nbsp;s (296&nbsp;MB) | **1.73&nbsp;s** (101&nbsp;MB) | 5.22&nbsp;s (311&nbsp;MB) | **1.88&nbsp;s** (100&nbsp;MB) |
| BGZF FASTQ | 35.1&nbsp;s | 23.4&nbsp;s | 6.75&nbsp;s (229&nbsp;MB) | **4.51&nbsp;s** (65&nbsp;MB) | 1.15&nbsp;s (228&nbsp;MB) | **0.82&nbsp;s** (185&nbsp;MB) | 0.56&nbsp;s (233&nbsp;MB) | **0.42&nbsp;s** (198&nbsp;MB) |
| SAM        | 47.1&nbsp;s | 19.6&nbsp;s | aborts* | **3.50&nbsp;s** (30&nbsp;MB) | aborts* | **0.60&nbsp;s** (94&nbsp;MB) | aborts* | **0.30&nbsp;s** (107&nbsp;MB) |
| BAM        | 38.6&nbsp;s | 17.2&nbsp;s | 8.04&nbsp;s (229&nbsp;MB) | **5.17&nbsp;s** (17&nbsp;MB) | 1.33&nbsp;s (215&nbsp;MB) | **0.98&nbsp;s** (129&nbsp;MB) | 1.48&nbsp;s (221&nbsp;MB) | **0.84&nbsp;s** (154&nbsp;MB) |

\* falco 2.0.3 aborts (exit 134) on this samtools-generated SAM file, which
samtools reads completely; falcoAmadeus, falco 1.3.2 and FastQC process all
reads.

**Paired-end, 10 million HG002 pairs of 2 x 126 bp (SRR12898355)**,
falcoAmadeus in its paired mode (`--paired-with`), the others on the two
files:

| Input | FastQC `-t 1` | falco 1.3.2 | falco 2.0.3 `-t 1` | falcoAmadeus `-t 1` | falco 2.0.3 `-t 8` | falcoAmadeus `-t 8` | falco 2.0.3 `-t 32` | falcoAmadeus `-t 32` |
| ----- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| FASTQ      | 49.2&nbsp;s | 21.7&nbsp;s | 6.55&nbsp;s (632&nbsp;MB) | **4.25&nbsp;s** (267&nbsp;MB) | 1.08&nbsp;s (566&nbsp;MB) | **0.68&nbsp;s** (314&nbsp;MB) | 0.44&nbsp;s (708&nbsp;MB) | **0.36&nbsp;s** (438&nbsp;MB) |
| gzip FASTQ | 62.1&nbsp;s | 32.7&nbsp;s | 18.3&nbsp;s (640&nbsp;MB) | **8.91&nbsp;s** (307&nbsp;MB) | 7.21&nbsp;s (572&nbsp;MB) | **3.16&nbsp;s** (464&nbsp;MB) | 6.83&nbsp;s (716&nbsp;MB) | **3.44&nbsp;s** (469&nbsp;MB) |

falcoAmadeus is faster than falco 2.0.3 with less memory in every cell, and
3 to 6 times faster than falco 1.3.2 on a single thread. The measurements,
including CPU time, cold page cache and runs pinned to 8 hardware threads,
are in [`paper/data`](paper/data/), and the figures are made by
[`paper/scripts/make_figures.py`](paper/scripts/make_figures.py).

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

The conda build links htslib, libdeflate and ISA-L, so it ships with full
performance on all input formats. The recipe is kept in
[`recipe/`](recipe/).

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

All of these are optional, detected at configure time, and fall back
gracefully when absent -- but the first two are strongly recommended, since
the fast parallel BAM and gzip paths depend on them:

- [htslib](https://github.com/samtools/htslib) is required to process BAM
  files (without it BAM is treated as an unrecognized format). Enable with
  `--enable-hts`.
- [libdeflate](https://github.com/ebiggers/libdeflate) powers parallel BGZF
  block decompression (BAM and bgzip-compressed FASTQ). Detected
  automatically when the header and library are installed.
- [ISA-L](https://github.com/intel/isa-l) roughly doubles plain-gzip
  decompression speed. Detected automatically.
- [rapidgzip](https://github.com/mxmlnkn/rapidgzip) (header-only) can
  decompress even a single plain-gzip stream in parallel. Build with
  `--with-rapidgzip=<checkout>/librapidarchive/src` and enable per run with
  `FA_RAPIDGZIP=1`; worthwhile for very large plain-gzip files on machines
  with spare cores and memory.

A typical full-performance build:

```
$ ./configure CXXFLAGS="-O3 -Wall" --enable-hts
```

with htslib, libdeflate and ISA-L installed where the compiler can find them
(add their locations to `CPPFLAGS`/`LDFLAGS` if they live in a non-default
prefix, e.g. a conda environment).

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

with copies allowed to differ in two bases, and an 8 nt UMI at the start of
read 2 that may carry one error:

```
$ falcoAmadeus R1.fastq.gz --paired-with R2.fastq.gz --dup-mismatches 2 \
    --umi-r2 1:8 --umi-mismatches 1 -o out/
```

**Assay-aware grading** (here for a bisulfite / EM-seq library):

```
$ falcoAmadeus reads.fastq.gz --assay wgbs -o out/
```

### Options

FalcoAmadeus-specific options:

| Option | Description |
| --- | --- |
| `--paired-with R2` | Analyse the R2 mate together with the R1 positional input; write one combined report and `paired_duplication.txt`. Rejects mates whose read identifiers disagree. |
| `--allow-unmatched-pairs` | With `--paired-with`, proceed even when R1/R2 read identifiers disagree. |
| `--dup-mismatches N` | Also report a duplicate rate in which copies may differ in up to N bases (at most 8): over both mates with `--paired-with` (`paired_duplication.txt`), per read otherwise (`aware_duplication.txt`). |
| `--umi-r1 START:LEN`, `--umi-r2 START:LEN` | The UMI lies in read 1 or read 2 at a 1-based position; pairs are copies only when their UMIs agree, and the UMI bases are left out of the sequence comparison. |
| `--umi-in-name` | The UMI is the last `:` or `_` separated field of the read 1 name (bcl-convert, UMI-tools extract, fastp). |
| `--optical` | Count optical duplicates (Picard's rule) within the groups of copies; reported in the aware duplication section. Needs Illumina read names (tile:x:y). |
| `--optical-distance N` | Pixel distance for `--optical` (default 100; Picard advises 2500 for patterned flow cells). |
| `--umi-mismatches N` | Bases in which the UMIs of two copies may differ (default 0). |
| `--acceptance FILE` | Grade the sample against acceptance rules; writes `acceptance.txt` and an `acceptance` section of the JSON summary. |
| `--acceptance-exit` | With `--acceptance`, exit with status 2 when a sample fails its rules. |
| `--screen DB` | Contamination screen against a database from `--screen-build`; writes `screen.txt`. |
| `--screen-reads N` | Reads per file for `--screen` (default 1000000; 0: all). |
| `--screen-build DB` with `--screen-refs name=file,...` | Build a screen database from FASTA files (at most 32 references), then exit. |
| `--trim-preview` | Write `trimming_preview.txt`: bases and reads removed by 3' quality trimming at 10, 20 and 30 (as cutadapt `-q`), reads left shorter than 20, 36 and 50 bases, and an adapter-trimming estimate. FASTQ files only. |
| `--trim-preview-reads N` | Reads per file for `--trim-preview` (default 1000000; 0: all). |
| `--insert-size-pairs N` | With `--paired-with`, the number of pairs (from the start) whose mate overlap gives `paired_insert_size.txt` (default 1000000; 0: all). |
| `--complexity` | Single-end: count every read's duplication key exactly and write `complexity.txt` with the expected distinct reads at 2x to 100x the sequenced reads (paired input reports it in `paired_duplication.txt` always). |
| `--paired-data mates\|combined\|both` | With `--paired-with`, which `fastqc_data.txt` files to emit for MultiQC: separate mates (default), a single pooled sample, or both. |
| `--assay TYPE` | Assay-aware pass/warn/fail grading (`rnaseq`, `amplicon`, `wgbs`/`emseq`/`rrbs`, `chipseq`, `atac`, `smallrna`, `wgs`, …; default `generic` is byte-identical to falco). |
| `--conversion-residual F` | For bisulfite/EM assays, the unconverted fraction outside CpG (0–1; C among C+T followed by a base other than G in C→T reads, G among G+A preceded by a base other than C in G→A reads; each read's strand from its composition) above which a possible incomplete-conversion warning is emitted (default 0.2). |
| `--legacy-n-content-grade` | Grade Per base N content as falco 1.x does (never warns or fails), for output byte-identical to falco 1.3.2. |
| `--name NAME` | Sample name for the report `Filename` / MultiQC, instead of the input basename. |
| `-t, --threads N` | Parse with N threads; works for FASTQ, gzipped FASTQ, SAM and BAM (byte-identical output at every thread count). |

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
