# Statistics beyond FastQC

The FastQC report describes single reads. The statistics below describe fragments and the library
as a whole. They are computed from the raw reads in the run that writes the FastQC report, and
where one ports a published method, it gives the same values as the original program on the same
reads. The measurements behind the statements on this page are in
[`paper/data`](https://github.com/joachimwolff/falcoAmadeus/tree/master/paper/data).

| Statistic | FastQC 0.12.1 | falco 1.3.2 | falco 2.0.3 | falcoAmadeus | Method ported or compared with |
|---|---|---|---|---|---|
| Duplicate rate of read pairs | no | no | no | exact | Picard MarkDuplicates (compared) |
| Duplicates differing by sequencing errors | no | no | no | up to 8 mismatches | Picard EstimateLibraryComplexity, clumpify (compared) |
| UMI-aware duplicates | no | no | no | in read or read name | UMI-tools, fgbio |
| Optical duplicates | no | no | no | yes | Picard OpticalDuplicateFinder (ported) |
| Library complexity at deeper sequencing | no | no | histogram for preseq | expected distinct reads or pairs | preseq lc_extrap (ported) |
| Insert size, adapter read-through, adapter dimers | no | no | no | from the mate overlap | fastp overlap analysis (ported) |
| Contamination by other genomes | no | no | no | k-mer screen | FastQ Screen, Kraken 2, Mash |
| Quality trimming preview | no | no | no | at Phred 10, 20, 30 | cutadapt quality trimming (ported) |
| UMI diversity, families and quality | no | no | no | yes | UMI-tools |
| Assay-specific grading | no | bisulfite | bisulfite | eleven assay types | |
| Acceptance rules per sample | no | no | no | yes | |
| JSON summary, MultiQC general statistics | no | no | no | yes | MultiQC |

## Paired duplicates

FastQC and falco read one file at a time, so a fragment that shares only one end with another
counts as a duplicate. falcoAmadeus keys each pair by both mates and counts the pairs that repeat
an earlier pair exactly, over all reads and independent of read order. On a simulated library with
known PCR duplication, this rate came closest to the true rate among seven programs.

## Duplicates with sequencing errors and UMIs

A copy with a sequencing error has a different key. With `--dup-mismatches N`, pairs whose mates
differ in at most N bases are grouped, found through blocks shared by both mates and confirmed base
by base; UMI options keep pairs with different UMIs apart. On public RNA-seq read pairs, two
mismatches covered 98 % of the duplicates that Picard MarkDuplicates marked by alignment position,
and 99 % of the counted pairs were Picard duplicates. On whole-genome sequencing, where many
duplicates differ in more bases, the coverage was about half.

## Optical duplicates

Within each group of copies, copies on the same tile within the pixel distance are optical
duplicates, with Picard's rule (one kept per cluster of close copies, the two mate orders apart,
coordinates in 16 bits). Where the groups equal Picard's, the counts equal Picard's.

## Library complexity

From the number of distinct pairs seen once, twice and so on, falcoAmadeus extrapolates the
distinct pairs expected at 2 to 100 times the sequenced pairs, as `preseq lc_extrap` does. The
curves equalled those of preseq 3.2.0 on the same counts.

## Insert size and adapter read-through

For the first million pairs, read 1 is aligned with the reverse complement of read 2 as fastp does.
The histogram equalled that of fastp 0.24.3 on the same pairs. The report adds the share of pairs
whose insert is shorter than a read, the adapter each read continues into and the adapter dimers.

## Contamination screen

Canonical 31-mers are selected by their hash value (FracMinHash) from the reads and from a panel of
reference genomes, and low-complexity k-mers are skipped. The report gives the share of read k-mers
found in each reference, found only there, and found in none.

## Trimming preview

The bases and reads that 3' quality trimming at Phred 10, 20 and 30 would remove, with cutadapt's
algorithm, and the reads left shorter than 20, 36 and 50 bases. The values equalled those of
cutadapt on public runs of eleven assays.

## UMI report

UMIs with N, homopolymer UMIs, distinct UMIs and the share of possible sequences in use, the most
frequent UMI and the base composition per position; the families of copies with their sizes; and
the base quality of UMIs in reads.

## Assay profiles

The generic FastQC limits flag the expected composition of many assays. The profiles change the
grades of the modules for eleven assay types. On 120 public runs of ten assays, they lowered the
flags on unmodified libraries for 16S, small RNA, RNA-seq and EM-seq, and every injected adapter,
quality, GC, duplicate and N defect was still detected. The profiles of 16S and small RNA suppress
the composition module and miss a composition bias there.
