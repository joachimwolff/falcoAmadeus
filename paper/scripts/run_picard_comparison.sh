#!/usr/bin/env bash
# Paired duplicate rates of falcoAmadeus and Picard MarkDuplicates on the
# same read pairs: align, mark duplicates, compare (compare_picard.py).
#
# Usage: run_picard_comparison.sh <R1.fastq> <R2.fastq> <bwa-mem2 index>
#                                 <falcoAmadeus> <out_dir> [threads]
# Needs bwa-mem2, samtools and gatk (MarkDuplicates) on PATH and a Python
# with numpy and pysam (PYTHON, default python3).
set -euo pipefail
R1=$1; R2=$2; REF=$3; FA=$4; OUT=$5; T=${6:-16}
HERE=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$OUT"

# bwa-mem2 requires equal mate names: drop the .1/.2 mate suffix of
# fastq-dump --readids names (SRR.N.1 / SRR.N.2 -> SRR.N)
for m in 1 2; do
  in=$R1; [ "$m" = 2 ] && in=$R2
  awk 'NR % 4 == 1 { sub(/\.[12]$/, "", $1); print $1; next } { print }' \
    "$in" > "$OUT/renamed_R$m.fastq"
done

bwa-mem2 mem -t "$T" -R '@RG\tID:pairs\tSM:pairs\tPL:ILLUMINA' "$REF" \
  "$OUT/renamed_R1.fastq" "$OUT/renamed_R2.fastq" 2> "$OUT/bwa.log" |
  samtools sort -@ 8 -m 2G -o "$OUT/aligned.bam" -
# duplicate-set tags (DI) only annotate; the duplicate flags are the same
gatk MarkDuplicates -I "$OUT/aligned.bam" -O "$OUT/marked.bam" \
  -M "$OUT/markduplicates_metrics.txt" --TAG_DUPLICATE_SET_MEMBERS true \
  > "$OUT/markduplicates.log" 2>&1
rm -f "$OUT/renamed_R1.fastq" "$OUT/renamed_R2.fastq" "$OUT/aligned.bam"

"${PYTHON:-python3}" -I "$HERE/compare_picard.py" "$OUT/marked.bam" "$R1" "$R2" "$OUT"
"${PYTHON:-python3}" -I "$HERE/picard_set_mismatches.py" "$OUT/marked.bam" \
  "$R1" "$R2" "$OUT/picard_set_mismatches.tsv"

# falcoAmadeus on all pairs and on the compared pairs only
for set in all compared; do
  a=$R1; b=$R2
  [ "$set" = compared ] && a=$OUT/compared_R1.fastq && b=$OUT/compared_R2.fastq
  mkdir -p "$OUT/fa_$set"
  "$FA" -t "$T" -q -o "$OUT/fa_$set" "$a" --paired-with "$b" > /dev/null
done
