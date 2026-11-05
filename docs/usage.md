# Usage

## Single-end input

```bash
falcoAmadeus -t 8 -o out/ reads.fastq.gz
```

The input can be FASTQ, gzip or BGZF FASTQ, SAM or BAM; `-` reads FASTQ from standard input
(`-f fastq` then gives the format). The run writes:

| File | Content |
|---|---|
| `fastqc_data.txt` | The FastQC report as text, read by MultiQC. |
| `summary.txt` | PASS, WARN or FAIL for each module. |
| `falcoAmadeus_report.html` | The interactive HTML report. |
| `falcoAmadeus_summary.json` | All metrics, including those of the optional analyses. |
| `falcoAmadeus_mqc.json` | General-statistics columns for MultiQC. |

## Paired-end input

```bash
falcoAmadeus -t 8 -o out/ R1.fastq.gz --paired-with R2.fastq.gz
```

The mates are analysed together. The run writes one HTML report with tabs for the combined view
(read 1 and read 2 as one data set) and each mate, the `fastqc_data.txt` of each mate for MultiQC
(`--paired-data` chooses mates, combined or both), and:

- `paired_duplication.txt`: the exact duplicate rate of read pairs, the per-mate rates, the FastQC
  estimates for comparison and the library complexity curve.
- `paired_insert_size.txt`: insert-size histogram, adapter read-through and adapter dimers from the
  first million pairs (`--insert-size-pairs`).

The read identifiers of the two files are compared, and files whose mates disagree are rejected
unless `--allow-unmatched-pairs` is given.

## Duplicates with sequencing errors and UMIs

```bash
falcoAmadeus -o out/ R1.fastq.gz --paired-with R2.fastq.gz \
    --dup-mismatches 2 --umi-r2 1:8 --umi-mismatches 1
```

Copies may differ in up to `--dup-mismatches` bases; with a UMI option (`--umi-r1`, `--umi-r2`
`START:LEN`, or `--umi-in-name`), pairs are copies only when their UMIs agree. A UMI option also
writes `umi_qc.txt`. For single-end input, the rate is written to `aware_duplication.txt`.

## Optical duplicates

```bash
falcoAmadeus -o out/ R1.fastq.gz --paired-with R2.fastq.gz --optical --optical-distance 2500
```

Counts the copies of one cluster read twice, with Picard's rule, from the tile and coordinates in
Illumina read names. Use 100 pixels for non-patterned and 2500 for patterned flow cells.

## Contamination screen

Build a database once from FASTA files, then screen runs against it:

```bash
falcoAmadeus --screen-build panel.db \
    --screen-refs human=GRCh38.fa,mouse=GRCm39.fa,phix=phix.fa,ecoli=ecoli.fa
falcoAmadeus -o out/ reads.fastq.gz --screen panel.db
```

`screen.txt` gives, for each reference, the share of read k-mers found in it, found only in it, and
found in none of the references.

## Trimming preview and complexity

```bash
falcoAmadeus -o out/ reads.fastq.gz --trim-preview --complexity
```

`trimming_preview.txt` gives the bases and reads that 3' quality trimming at Phred 10, 20 and 30
would remove, as cutadapt `-q` does, and an adapter-trimming estimate. `--complexity` writes
`complexity.txt` for single-end input; paired input always reports complexity.

## Assay profiles and acceptance rules

```bash
falcoAmadeus -o out/ reads.fastq.gz --assay wgbs --acceptance rules.tsv
```

`--assay` grades the modules with the expectations of the assay: `wgs`, `exome`, `amplicon`, `16s`,
`rnaseq`, `smallrna`, `atac`, `chipseq`, `wgbs`, `emseq` or `rrbs` (default `generic`, the FastQC
limits). The reported values do not change, only the grades. An acceptance file holds one rule per
line, for example:

```text
q30_bases_pct                          >=  85  75
screen.references.human.found_pct      >=  90  80
```

The metric is a MultiQC column or a path into the JSON summary, followed by the WARN and FAIL
limits. `acceptance.txt` lists each rule and an overall PASS, WARN or FAIL;
`--acceptance-exit` makes a failed sample end the run with status 2.
