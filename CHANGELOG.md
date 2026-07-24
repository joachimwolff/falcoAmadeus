# Changelog

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
