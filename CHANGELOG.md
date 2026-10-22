# Changelog

## FalcoAmadeus 4 (unreleased)

### Changed

- **Paired-end duplication is counted exactly.** `paired_duplication.txt`
  previously reported the FastQC estimator (which tracks the first 100,000
  distinct sequences and extrapolates) as the paired rate, so the result
  depended on where duplicates sat in the file: 45% true duplication read as
  0% when the copies came late and 47.4% when early. The paired, Read1-only
  and Read2-only rates are now counted exactly over all reads from 64-bit
  fingerprints and are independent of read order. The FastQC-style estimates
  remain as separately labelled rows for comparison with FastQC/falco. The
  row names changed (`Paired-end exact`, `Paired-end FastQC-style`, ...).
- **Idle threads sleep instead of spinning** (`OMP_WAIT_POLICY=passive`
  unless the user sets it): 10-35% less CPU time at unchanged wall time; on
  an 8-CPU machine `-t 8` now uses 15% more CPU than `-t 1` (was 60%).

### Fixed

- **Mates out of order are rejected.** `--paired-with` accepted R1/R2 files
  whose read identifiers did not correspond and reported a duplicate rate
  for meaningless pairs. Mismatched identifiers now stop the run with an
  explanatory error; `--allow-unmatched-pairs` overrides.

### Performance

Measured on a 32-core Threadripper 7970X (whole machine, and pinned to 4
cores / 8 hardware threads), cold and warm page cache, 1 to 32 threads, on
public data (GSM4212622 10M reads; SRR14777166 29.9M reads as FASTQ, gzip,
BGZF, SAM and BAM). falcoAmadeus is now faster than falco 2.0.3 and uses
less memory in all 96 time and memory comparisons (falco 2.0.3 aborts on
the SAM file); before these changes it lost 14.

- **No cache-line sharing between parsing threads.** The per-range GC
  record vectors and per-thread tile record vectors were adjacent objects
  appended once per read by neighbouring threads. Each range now appends on
  the parsing thread's stack. 29.9M reads, `-t 32`: 0.43 s (falco 2.0.3)
  vs 0.53 s before, 0.33 s now.
- **Cold-cache reads.** Uncompressed FASTQ that is not in the page cache is
  read with sequential `pread` into the parse pipeline instead of being
  memory-mapped (kernel readahead stays effective); the BGZF block index is
  walked on a background thread instead of reading the whole compressed
  file before the first chunk. BGZF `-t 8`, cold: 1.23 s to 0.90 s.
- **Parallel BGZF FASTQ.** Block groups are read and inflated by the parse
  team itself, overlapping the parse of the previous chunk.
- **Lower fixed costs:** page-cache residency sampled at 1024 pages instead
  of a full `mincore` walk (~70 ms on a 4 GB file); FastQC's duplication
  extrapolation runs 8 duplication levels per AVX2 pass, bit-identical
  (0.19 s to 0.075 s); GC models built lazily per read length (~9 ms); HTML
  placeholders replaced without `std::regex` (~12 ms); no passive-wait
  re-exec at `-t 1`; range boundaries found by the parsing threads.
- **Memory:** per-thread duplication counts in 32-bit cells (folded into
  the map before they could overflow), per-position quality tables in
  calloc'd zero pages, 64 MB memory-mapped slices, and at most 8 threads
  with 32 MB chunks for plain gzip, whose single serial inflater bounds the
  throughput. Peak RSS at `-t 32`: FASTQ 179 MB (falco 2.0.3: 326 MB), gzip
  111 MB (308 MB), BGZF 215 MB (236 MB), BAM 182 MB (222 MB).
- Paired mode reads both mates concurrently in large ISA-L/zlib blocks and
  streams the pooled combined view instead of writing a concatenated
  temporary copy of both inputs (tiled input and builds without ISA-L keep
  the temporary copy). 2M simulated pairs, `-t 8`: uncompressed 2.64 s to
  0.93 s, gzip 5.5 s to 2.8 s.

### Tests

- `test/paired/`: a deterministic paired-end simulator with known duplicate
  fractions and a test script (in CI) covering late/early/shuffled copies,
  shared-R1-only copies, mate reordering and thread invariance.
- The upstream parity gate now compares against falco 1.3.2 (identical
  output to 1.3.1 on the bundled tests).

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
  (~36x falco 1.3.1), gzipped FASTQ 40&nbsp;s (~10x). The unreleased falco2
  master of 2026-08-05 could not be built with the conda clang/libc++
  toolchain on macOS (C++23 `std::views::enumerate`); falco 2.0.2 and later
  add macOS build support.

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
