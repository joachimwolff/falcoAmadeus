# Changelog

## FalcoAmadeus 3

A portability and memory release, verified end to end on Apple Silicon
(M3 Pro, macOS 26.6) with byte-identical output to the Linux build.

### Fixed

- **macOS builds compile again.** Two `madvise(..., MADV_HUGEPAGE)` calls
  (Linux-only optimization hints, introduced with the v2 BAM work) broke the
  macOS arm64 build at compile time; both are now compiled conditionally.
- **BAM memory on macOS.** The BGZF block-index walk memory-mapped the whole
  compressed file and relied on `MADV_DONTNEED` drop-behind, which Darwin
  treats as a lazy hint: a 4&nbsp;GB BAM peaked at 3.1 to 3.8&nbsp;GB resident
  on macOS. The walk now streams through a bounded 8&nbsp;MB `pread` window:
  ~150&nbsp;MB peak on macOS (unchanged on Linux), and BAM runs about 30%
  faster on Apple Silicon as a side effect.

### Added

- **macOS CI.** The push/PR pipeline now includes a macos job building with
  the full optional dependency set (htslib, libdeflate, ISA-L), so
  Darwin-only breakage is caught before release.

### Notes

- On Apple Silicon (163M reads, M3 Pro): BAM 9.7&nbsp;s at `-t 12`
  (~36x falco 1.3.1), gzipped FASTQ 40&nbsp;s (~10x). falco2 cannot be built
  on macOS at all (its source needs C++23 `std::views::enumerate`, which
  libc++ does not implement).

## FalcoAmadeus 2

A performance release: every input format now parses in parallel, and
single-threaded runs are faster too. Output remains **byte-identical** to
upstream falco at every thread count; each change in this release was gated on
the upstream md5 test suite, on diffs against the previous binary across all
formats (including high-duplication and reverse-strand BAM inputs), and on
`-t 1 == -t 8 == -t 16` output equality.

### Added

- **Parallel parsing for every input format** (`-t N`): uncompressed FASTQ and
  SAM are memory-mapped and parsed in place across record-aligned byte ranges;
  BAM decompresses and parses in overlapped, double-buffered chunks (BGZF
  blocks inflate in parallel via libdeflate); plain gzip inflates with ISA-L
  into bounded chunks parsed in parallel. Single-core runs use the same
  drivers, so `-t 1` is faster as well.
- **Standard input** (`-`): all formats can be read from a pipe
  (`--format` required; the per-tile module is unavailable on pipes).
- **ISA-L igzip** support for plain-gzip inflate (autodetected at configure
  time, optional).
- **rapidgzip backend** (optional, `--with-rapidgzip=DIR` at configure time,
  enabled per run with `FA_RAPIDGZIP=1`): parallel decompression of large
  plain-gzip files on machines with spare cores and memory.
- README benchmark section with the 163M-read three-tool comparison
  (falco 1.3.1, falco2, FalcoAmadeus) and the accompanying figure.

### Performance

Measured on 163,264,569 reads (AMD Ryzen 9 7950X, warm page cache),
FalcoAmadeus at `-t 16`: FASTQ 3.6&nbsp;s, gzipped FASTQ 17.3&nbsp;s, SAM
5.7&nbsp;s, BAM 6.1&nbsp;s. That is 12x (gzip) to 53x (FASTQ) faster than
single-threaded falco 1.3.1 and faster than falco2 in every format at both
`-t 8` and `-t 16`, at lower peak memory (133 to 239&nbsp;MB vs falco2's 211
to 308&nbsp;MB). The main mechanisms:

- SIMD statistics kernels (AVX-512, SSSE3 and BMI2 on x86, NEON on ARM,
  scalar fallback everywhere): the sequence pass, the quality pass with
  cache-resident accumulator tiles, chain-free adapter window probes, and
  BAM packed-base decode including reverse-complement.
- Order-exact parallel bookkeeping: duplication's first-100k-unique
  semantics, the floating-point GC model and read subsampling reproduce the
  serial trajectory exactly, which is what keeps output byte-identical.
- Duplication tracker: open-addressing map keyed directly off the sequence
  buffer, hash-fingerprint slot tags with a flat shadow key store, lazy
  allocation in parallel workers, and prefetched split-phase lookups.
- Bounded memory: inputs stream through reused chunk buffers with
  drop-behind page management and cold-file readahead, so peak memory stays
  flat (roughly 35 to 250&nbsp;MB by format and thread count) regardless of
  input size.

### Fixed

- Default-adapter runs (no `adapter_list.txt` in the working directory,
  i.e. every typical run) silently used the slow string-search matcher
  instead of the optimized one; this upstream falco bug roughly doubled
  whole-run cost on all formats.
- `--version` printed `falco 1`; it now prints `falcoAmadeus 2`.
- The overrepresented-sequences listing is sorted deterministically
  (frequency, then sequence) regardless of hash-map iteration order.
- Paired-view HTML report no longer renders a duplicate top bar; the module
  sidebar is collapsible.
- `make distcheck`: `src/dup_hash_map.hpp` is shipped in release tarballs,
  committed test-suite artifacts were removed (with a dist-hook guard), and
  VPATH builds get a `test/` log directory.

### Changed

- Version 2: the `fastqc_data.txt` header line is now `##FalcoAmadeus 2` and
  the test-suite checksums were regenerated for it (outputs differ from
  version 1 only in that line).

## FalcoAmadeus 1

First release of **FalcoAmadeus**, a performance-focused and sequencing-aware
fork of [falco](https://github.com/smithlabcode/falco) (a fast reimplementation
of FastQC). A plain single-end run's `fastqc_data.txt` and `summary.txt` are
**byte-identical** to upstream falco (the upstream md5 test suite passes);
everything new below is **opt-in**.

### Added

- **Correct paired-end analysis** (`--paired-with R2`): a single combined HTML
  report with **Combined / Forward / Reverse** tabs (Combined pools R1+R2 as one
  dataset); the **corrected paired PCR-duplicate rate** — a pair is a duplicate
  only when *both* mates are — and how many R1/R2 reads **match** into pairs.
  Per-mate `fastqc_data.txt` for MultiQC plus a `paired_duplication.txt` summary.
- **`--paired-data mates|combined|both`**: choose whether MultiQC gets the two
  mates separately, a single pooled sample per pair, or both.
- **Assay-aware QC** (`--assay <type>`): per-library-type pass/warn/fail grading
  (`rnaseq`, `amplicon`, `wgbs`/`emseq`/`rrbs`, `chipseq`, `atac`, `smallrna`,
  `wgs`, …), so expected assay biology is not mis-flagged as failure.
- **`--conversion-residual`**: bisulfite incomplete-conversion warning threshold.
- **`--name NAME`**: explicit sample name for the report `Filename` / MultiQC
  (fixes the *forward*/*reverse* sample collapse for Galaxy paired collections).
- **Redesigned interactive HTML report** (`falcoAmadeus_report.html`):
  light/dark theme toggle (dark by default), selectable colourblind-safe colour
  palettes that also recolour the box plot and per-tile heatmap, collapsible
  modules with plain-language explanations and FastQC documentation links,
  readable unified hover tooltips, and cleaner plots (box-and-whisker per-base
  quality over pass/warn/fail zones, padded length-distribution histogram).
- **Multi-threaded gzip parsing** (`-t N`): background decompression overlapped
  with parsing, byte-identical to the single-threaded default.
- A **Galaxy tool wrapper** (maintained as `tools-iuc/tools/falcoAmadeus`).

### Performance

- Input-parser optimizations (`src/StreamReader.*`): roughly **1.2–1.4× faster**
  single-threaded on FASTQ/FASTQ.gz/SAM; peak memory unchanged (input is still
  streamed).
- Multi-threaded gzip parsing (`-t N`) reaches **~2.5× over single-thread** on
  gzipped input (byte-identical), i.e. roughly **3× faster than upstream falco**
  on gzipped data. Uncompressed input is unaffected.

### Changed

- Executable renamed to **`falcoAmadeus`**; HTML report named
  `falcoAmadeus_report.html`; the `fastqc_data.txt` identifier is
  `##FalcoAmadeus` (still ingested by MultiQC as a FastQC report).
- Removed the FastQC compatibility dummy options that falco silently ignored
  (`--casava`, `--nano`, `--nofilter`, `--extract`, `--java`, `--noextract`,
  `--min_length`, `--kmers`, `--dir`).

### Fixed

- MultiQC can now ingest `--assay` reports (assay comment lines are emitted
  between modules instead of inside them, where MultiQC's parser choked).
- CI: `autogen.sh` no longer requires the checkout directory to be named
  literally `falco`, so it works for the `falcoAmadeus` fork.

### Compatibility

- Default single-end `fastqc_data.txt` and `summary.txt` remain byte-identical
  to upstream falco.
