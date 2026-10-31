# UMI quality on a public UMI run

SRR27937408 (PRJNA1075320; Lexogen QuantSeq 3' mRNA-Seq FWD with the UMI
Second Strand Synthesis Module, single-end 100 bases, 4,987,876 reads; ENA
MD5 2f841f16665789c62281eba11bef3b25): read 1 starts with a 6-nt UMI and the
spacer TATA.

`SRR27937408_mismatchesD.umi_qc.tsv` is umi_qc.txt of

    falcoAmadeus -t 16 SRR27937408.fastq.gz --umi-r1 1:6 --dup-mismatches D

for D = 0, 1, 2 (the File line shortened to the file name). For each D the
file equals `test/umi_qc/umi_oracle.py` run on the same reads and on the
groups of the same run (FA_AWARE_CLUSTERS), line for line.
