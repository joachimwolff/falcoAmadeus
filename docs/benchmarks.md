# Benchmarks

Measured on an AMD Ryzen Threadripper 7970X (32 cores, 64 threads) with a warm page cache, using
the whole machine; medians of two (single-end) or three (paired-end) runs. falco 2.0.3 was run
with its default duplication and with `--orig-dups`, and the faster of the two is shown. falco 1.3.2
and FastQC run single-threaded. Wall time, with peak memory in parentheses.

## Single-end

29.9 million reads of 1 to 76 bases (SRR14777166).

| Input | FastQC | falco 1.3.2 | falco 2.0.3 `-t 1` | falcoAmadeus `-t 1` | falco 2.0.3 `-t 8` | falcoAmadeus `-t 8` | falco 2.0.3 `-t 32` | falcoAmadeus `-t 32` |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| FASTQ | 32.6 s | 20.4 s | 5.70 s (287 MB) | **3.59 s** (31 MB) | 1.04 s (299 MB) | **0.57 s** (95 MB) | 0.49 s (322 MB) | **0.23 s** (103 MB) |
| gzip FASTQ | 34.6 s | 23.2 s | 9.65 s (298 MB) | **4.35 s** (63 MB) | 5.28 s (296 MB) | **1.73 s** (101 MB) | 5.22 s (311 MB) | **1.88 s** (100 MB) |
| BGZF FASTQ | 35.1 s | 23.4 s | 6.75 s (229 MB) | **4.51 s** (65 MB) | 1.15 s (228 MB) | **0.82 s** (185 MB) | 0.56 s (233 MB) | **0.42 s** (198 MB) |
| SAM | 47.1 s | 19.6 s | aborts | **3.50 s** (30 MB) | aborts | **0.60 s** (94 MB) | aborts | **0.30 s** (107 MB) |
| BAM | 38.6 s | 17.2 s | 8.04 s (229 MB) | **5.17 s** (17 MB) | 1.33 s (215 MB) | **0.98 s** (129 MB) | 1.48 s (221 MB) | **0.84 s** (154 MB) |

falco 2.0.3 aborts (exit 134) on this samtools-generated SAM file, which samtools reads completely.

## Paired-end

10 million HG002 pairs of 2 x 126 bases (SRR12898355); falcoAmadeus in its paired mode, the other
programs on the two files.

| Input | FastQC | falco 1.3.2 | falco 2.0.3 `-t 1` | falcoAmadeus `-t 1` | falco 2.0.3 `-t 8` | falcoAmadeus `-t 8` | falco 2.0.3 `-t 32` | falcoAmadeus `-t 32` |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| FASTQ | 49.2 s | 21.7 s | 6.55 s (632 MB) | **4.25 s** (267 MB) | 1.08 s (566 MB) | **0.68 s** (314 MB) | 0.44 s (708 MB) | **0.36 s** (438 MB) |
| gzip FASTQ | 62.1 s | 32.7 s | 18.3 s (640 MB) | **8.91 s** (307 MB) | 7.21 s (572 MB) | **3.16 s** (464 MB) | 6.83 s (716 MB) | **3.44 s** (469 MB) |

falcoAmadeus was faster than falco 2.0.3 with less memory in every configuration, also with a cold
page cache and when restricted to 8 hardware threads, and 3 to 6 times faster than falco 1.3.2 on a
single thread. The measurements, including CPU time, are in
[`paper/data`](https://github.com/joachimwolff/falcoAmadeus/tree/master/paper/data), and the
figures are made by
[`paper/scripts/make_figures.py`](https://github.com/joachimwolff/falcoAmadeus/blob/master/paper/scripts/make_figures.py).

## How the speed is reached

- **Parallel parsing of every format.** Uncompressed FASTQ and SAM are memory-mapped and parsed in
  record-aligned ranges; BGZF blocks of BAM and bgzip FASTQ are decompressed in parallel with
  libdeflate; plain gzip is decompressed with ISA-L into chunks parsed in parallel while the next
  chunk is decompressed. Single-threaded runs use the same drivers.
- **Vector instructions** (AVX-512 on x86, NEON on ARM, scalar otherwise) for base composition,
  quality, adapter search and BAM base decoding.
- **Order-exact bookkeeping.** Everything that depends on read order in falco (the duplication
  estimate from the first 100,000 distinct sequences, the GC model, subsampling) reproduces the
  serial result, so the output is identical at every thread count.
- **Bounded memory.** Inputs stream through reused buffers, so peak memory does not grow with the
  size of the input.
