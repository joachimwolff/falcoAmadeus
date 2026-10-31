#!/usr/bin/env bash
# Library complexity and insert size of falcoAmadeus against the programs
# whose methods they port, on public read pairs: the yield curve against
# preseq 3.2.0 lc_extrap -Q on the same copy-number histogram, and the
# insert-size histogram against fastp 0.24.3 (one worker thread, so that
# fastp counts every pair) on all pairs.
#
# Usage: run_complexity_insert.sh <falcoAmadeus> <preseq> <fastp> <out_dir>
#                                 <name> <R1.fastq> <R2.fastq> [...]
# (name, R1, R2 repeated per data set). Appends to <out_dir>/summary.tsv
# and writes <out_dir>/<name>.yield.tsv (falcoAmadeus's complexity rows).
set -euo pipefail
FA=$1; PQ=$2; FP=$3; OUT=$4; shift 4
mkdir -p "$OUT"
S=$OUT/summary.tsv
[ -f "$S" ] || printf 'dataset\tpairs\tcheck\tvalues_compared\tequal\n' > "$S"
while [ $# -ge 3 ]; do
  name=$1; r1=$2; r2=$3; shift 3
  W=$OUT/work_$name; rm -rf "$W"; mkdir -p "$W"
  FA_COPY_HIST=$W/hist.txt FA_COMPLEXITY_CURVE=$W/curve.txt \
    "$FA" -t 16 -q --insert-size-pairs 0 -o "$W/fa" "$r1" \
    --paired-with "$r2" > "$W/fa.log" 2>&1
  pairs=$(awk -F'\t' '$1 == "Total read pairs" {print $2; exit}' \
    "$W/fa/paired_duplication.txt")
  "$PQ" lc_extrap -Q -H "$W/hist.txt" -o "$W/preseq.txt" 2> "$W/preseq.log"
  n=$(($(wc -l < "$W/preseq.txt") - 1))
  eq=$(cmp -s "$W/curve.txt" "$W/preseq.txt" && echo yes || echo no)
  printf '%s\t%s\tyield curve vs preseq lc_extrap -Q\t%s\t%s\n' \
    "$name" "$pairs" "$n" "$eq" >> "$S"
  awk '/^#Library complexity/ {f = 1} f && /^(#Sequenced|[0-9]+x\t)/' \
    "$W/fa/paired_duplication.txt" > "$OUT/$name.yield.tsv"
  "$FP" -w 1 -A -G -i "$r1" -I "$r2" -o /dev/null -O "$W/fp2.fq" \
    -j "$W/fp.json" -h /dev/null 2> "$W/fastp.log"
  rm -f "$W/fp2.fq"
  python3 -I -c 'import json,sys; j=json.load(open(sys.argv[1]))["insert_size"]
for v in j["histogram"] + [j["unknown"]]: print(v)' "$W/fp.json" > "$W/fp.txt"
  awk -F'\t' '$1 ~ /^[0-9]+$/ && NF == 2 {h[$1] = $2} $1 ~ /^>=512/ {h[512] = $2}
    END {for (i = 0; i <= 512; i++) print (i in h) ? h[i] : 0}' \
    "$W/fa/paired_insert_size.txt" > "$W/fa_hist.txt"
  eq=$(cmp -s "$W/fa_hist.txt" "$W/fp.txt" && echo yes || echo no)
  printf '%s\t%s\tinsert-size histogram vs fastp\t%s\t%s\n' \
    "$name" "$pairs" "$(wc -l < "$W/fp.txt")" "$eq" >> "$S"
  cp "$W/fa/paired_insert_size.txt" "$OUT/$name.insert_size.txt"
done
