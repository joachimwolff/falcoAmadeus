# Trimming preview against cutadapt

`cutadapt_comparison.tsv` is the output of
`paper/scripts/trim_preview_validation/compare_cutadapt.py` (falcoAmadeus 4,
cutadapt 5.2) on the first 200,000 reads of 11 public runs:

| file | run | assay |
|---|---|---|
| rnaseq_SRR41451237 | SRR41451237 | RNA-seq |
| smallrna_SRR41413611 | SRR41413611 | small RNA |
| atac_SRR36625746 | SRR36625746 | ATAC-seq |
| chipseq_SRR40896959 | SRR40896959 | ChIP-seq |
| wgbs_SRR39134064 | SRR39134064 | WGBS |
| emseq_SRR13076814 | SRR13076814 | EM-seq |
| rrbs_SRR41080386 | SRR41080386 | RRBS |
| amplicon16s_SRR1619438 | SRR1619438 | 16S amplicon |
| exome_SRR41127230 | SRR41127230 | exome |
| wgs_SRR14777166 | SRR14777166 | WGS |
| pe_SRR12898355_R1 | SRR12898355 (read 1) | WGS, paired |

Rows with cutoff 10, 20 and 30 compare the preview's quality rows with the
values counted from `cutadapt -q CUTOFF` (column `equal`: all 165 True).
Rows `adapter+20` compare the adapter estimate with `cutadapt -a ADAPTER
-q 20`, ADAPTER being the 12 bases the preview chose (column `estimate`).
