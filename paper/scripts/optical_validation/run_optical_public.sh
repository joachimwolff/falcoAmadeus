#!/usr/bin/env bash
# Optical duplicates of falcoAmadeus (--optical, alignment-free) against
# Picard MarkDuplicates (GATK 4.6.2.0) on public runs with Illumina read
# names: the first N pairs are aligned with bwa-mem2 (hg19), marked with
# Picard at 100 and 2500 pixels, and counted by falcoAmadeus with 0, 1 and
# 2 mismatches; FA_AWARE_CLUSTERS keeps falcoAmadeus's groups (100 pixels)
# for compare_on_picard_pairs.py.
# usage: run_optical_public.sh FALCOAMADEUS OUTDIR NAME R1.fastq.gz R2.fastq.gz PAIRS
# needs bwa-mem2, samtools and gatk on PATH; REF (default ~/transvar_test/hg19.fa)
set -eu
FA=$1 OUT=$2 NAME=$3 R1=$4 R2=$5 PAIRS=$6
REF=${REF:-$HOME/transvar_test/hg19.fa}
d=$OUT/$NAME
mkdir -p "$d"
zcat "$R1" | head -n $((PAIRS * 4)) > "$d/r1.fq"
zcat "$R2" | head -n $((PAIRS * 4)) > "$d/r2.fq"
bwa-mem2 mem -t 24 -R '@RG\tID:rg1\tSM:s1\tLB:lib1' "$REF" "$d/r1.fq" "$d/r2.fq" \
  2> "$d/bwa.log" | samtools sort -@8 -m 1500M -o "$d/aln.bam" - 2> "$d/sort.log"
samtools view -f 65 -F 2316 "$d/aln.bam" |
  awk '{print $1 "\t" $3 ":" $4 ":" $8 ":" int($2 / 16) % 2 "\t" $5}' > "$d/picard_pairs.tsv"
for D in 100 2500; do
  gatk --java-options -Xmx16g MarkDuplicates -I "$d/aln.bam" -O "$d/md$D.bam" \
    -M "$d/md$D.txt" --OPTICAL_DUPLICATE_PIXEL_DISTANCE $D \
    --VALIDATION_STRINGENCY SILENT > "$d/md$D.log" 2>&1
  rm -f "$d/md$D.bam"
  for M in 0 1 2; do
    rm -rf "$d/fa_${D}_$M"
    FA_AWARE_CLUSTERS="$d/groups_${D}_$M.bin" /usr/bin/time -f "%e s %M KB" \
      -o "$d/fa_${D}_$M.time" "$FA" -q -t 16 -o "$d/fa_${D}_$M" "$d/r1.fq" \
      --paired-with "$d/r2.fq" --optical --optical-distance $D --dup-mismatches $M
  done
done
