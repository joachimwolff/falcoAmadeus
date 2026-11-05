# Options

falcoAmadeus keeps the command line of falco; the options added by falcoAmadeus are listed below.
Run `falcoAmadeus --help` for the complete list.

## Input and performance

| Option | Meaning |
|---|---|
| `-t, --threads N` | Parse with N threads; FASTQ, gzip and BGZF FASTQ, SAM and BAM. The output is identical at every thread count. |
| `-f, --format FMT` | Force the input format (`fastq`, `fastq.gz`, `sam`, `bam`, ...); needed for standard input (`-`). |
| `--name NAME` | Sample name for the report and MultiQC, instead of the input file name. |

## Paired-end input

| Option | Meaning |
|---|---|
| `--paired-with R2` | Analyse read 2 together with the read 1 input: one combined report, `paired_duplication.txt` and `paired_insert_size.txt`. |
| `--allow-unmatched-pairs` | Proceed when the read identifiers of the two files disagree. |
| `--paired-data mates\|combined\|both` | Which `fastqc_data.txt` files to write for MultiQC: the two mates (default), the combined view, or both. |
| `--insert-size-pairs N` | Pairs from the start used for the insert-size report (default 1000000; 0: all). |

## Duplicates

| Option | Meaning |
|---|---|
| `--dup-mismatches N` | Also report a duplicate rate in which copies may differ in up to N bases (at most 8), over both mates for paired input. |
| `--umi-r1 START:LEN`, `--umi-r2 START:LEN` | The UMI lies in read 1 or read 2 at a 1-based position; copies must agree in their UMIs, and the UMI bases are left out of the sequence comparison. Also writes `umi_qc.txt`. |
| `--umi-in-name` | The UMI is the last `:` or `_` separated field of the read 1 name. |
| `--umi-mismatches N` | Bases in which the UMIs of two copies may differ (default 0). |
| `--optical` | Count optical duplicates with Picard's rule; needs Illumina read names with tile and coordinates. |
| `--optical-distance N` | Pixel distance for `--optical` (default 100; 2500 for patterned flow cells). |
| `--complexity` | Single-end: write `complexity.txt` with the expected distinct reads at 2 to 100 times the sequencing depth. |

## Contamination screen and trimming preview

| Option | Meaning |
|---|---|
| `--screen DB` | Screen the reads against a database built with `--screen-build`; writes `screen.txt`. |
| `--screen-reads N` | Reads per file screened (default 1000000; 0: all). |
| `--screen-build DB` | Build a screen database from `--screen-refs`, then exit. |
| `--screen-refs name=fasta,...` | FASTA files of the references (at most 32). |
| `--trim-preview` | Write `trimming_preview.txt`: 3' quality trimming at Phred 10, 20 and 30 as cutadapt `-q`, and an adapter-trimming estimate. FASTQ input only. |
| `--trim-preview-reads N` | Reads per file for the preview (default 1000000; 0: all). |

## Grading

| Option | Meaning |
|---|---|
| `--assay TYPE` | Assay profile: `generic` (default), `wgs`, `exome`, `amplicon`, `16s`, `rnaseq`, `smallrna`, `atac`, `chipseq`, `wgbs`, `emseq`, `rrbs`. |
| `--conversion-residual F` | For `wgbs`, `emseq` and `rrbs`, the unconverted fraction outside CpG above which a possible incomplete conversion is reported (default 0.2). |
| `--acceptance FILE` | Grade the sample against acceptance rules; writes `acceptance.txt`. |
| `--acceptance-exit` | With `--acceptance`, exit with status 2 when the sample fails. |
| `--legacy-n-content-grade` | Grade Per base N content as falco 1.x does, for output identical to falco 1.3.2. |

## Options kept from falco

`-o/--outdir`, `-c/--contaminants`, `-a/--adapters`, `-l/--limits`, `--nogroup`,
`-s/-subsample`, `-b/-bisulfite`, `-r/-reverse-complement`, the `-skip-*` options and the
`-D`/`-R`/`-S` file names work as in falco. The FastQC options that falco accepted and ignored
(`--casava`, `--nano`, `--extract`, `--java`, `--kmers`, `--dir`, `--min_length`) were removed.
