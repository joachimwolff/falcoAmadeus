# Optical duplicates against Picard on public runs

Two public paired runs whose FASTQ files keep the Illumina read names
(tile:x:y):

| dataset | source | instrument | pairs used |
|---|---|---|---|
| hiseq2500_HG002 | GIAB HG002 NIST_Illumina_2x250bps D1_S1_L001 (NCBI FTP over HTTPS) | HiSeq 2500 rapid, 2x250 (trimmed, 35-250), PCR-free | first 5,000,000 |
| novaseq_ERR5654681 | ENA ERR5654681 submitted files (PRJEB42052, A549 input) | NovaSeq 6000, 2x51 | first 10,000,000 |

`paper/scripts/optical_validation/run_optical_public.sh` aligns the pairs
with bwa-mem2 2.2.1 to hg19, marks duplicates with Picard MarkDuplicates
(GATK 4.6.2.0) at 100 and 2500 pixels, and runs falcoAmadeus `--optical`
with 0, 1 and 2 mismatches (`summary.tsv`; seconds and peak memory of
falcoAmadeus at -t 16 on an AMD Ryzen 9 7950X workstation).

Picard examines aligned pairs only (6,270,924 of the NovaSeq pairs, many of
whose inserts are shorter than the 30 bases bwa mem needs); falcoAmadeus
leaves out pairs that would not align (`falco_left_out`) and counts the
rest. `compare_on_picard_pairs.py` restricts falcoAmadeus's groups to the
pairs Picard examines, splits them by Picard's orientation and counts
optical duplicates by Picard's rule (`on_picard_pairs.tsv`). It also counts
close copies (100 pixels) within a falcoAmadeus group that Picard placed at
different coordinates, and how many of them have mapping quality 20 or more
in both pairs: on NovaSeq 22,926 of 46,836 such copies at 1 mismatch, 836
with mapping quality 20 or more, i.e. repeat reads aligned to different
copies, which Picard puts into different sets.
