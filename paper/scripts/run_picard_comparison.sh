#!/usr/bin/env bash
# Paired duplicate rates of falcoAmadeus and Picard MarkDuplicates on the
# same read pairs: align, mark duplicates, compare (compare_picard.py).
#
# Usage: run_picard_comparison.sh <R1.fastq> <R2.fastq> <index>
#                                 <falcoAmadeus> <out_dir> [threads]
# ALIGNER=bwa (default; index = bwa-mem2 index prefix, for DNA) or
# ALIGNER=star (index = STAR genome directory, for RNA-seq).
# Needs the aligner, samtools and gatk (MarkDuplicates) on PATH and a Python
# with numpy and pysam (PYTHON, default python3).
set -euo pipefail
R1=$1; R2=$2; REF=$3; FA=$4; OUT=$5; T=${6:-16}
HERE=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$OUT"

# the aligners require equal mate names: both mates of the n-th pair are
# named p.n (the comparison scripts read the pair index from it)
for m in 1 2; do
  in=$R1; [ "$m" = 2 ] && in=$R2
  awk 'NR % 4 == 1 { print "@p." (NR + 3) / 4; next } { print }' \
    "$in" > "$OUT/renamed_R$m.fastq"
done

if [ "${ALIGNER:-bwa}" = star ]; then
  STAR --runThreadN "$T" --genomeDir "$REF" \
    --readFilesIn "$OUT/renamed_R1.fastq" "$OUT/renamed_R2.fastq" \
    --outSAMtype BAM SortedByCoordinate --outSAMattrRGline ID:pairs SM:pairs \
    --outFileNamePrefix "$OUT/star_" > "$OUT/star.log" 2>&1
  mv "$OUT/star_Aligned.sortedByCoord.out.bam" "$OUT/aligned.bam"
else
  bwa-mem2 mem -t "$T" -R '@RG\tID:pairs\tSM:pairs\tPL:ILLUMINA' "$REF" \
    "$OUT/renamed_R1.fastq" "$OUT/renamed_R2.fastq" 2> "$OUT/bwa.log" |
    samtools sort -@ 8 -m 2G -o "$OUT/aligned.bam" -
fi
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
