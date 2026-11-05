#!/usr/bin/env bash
# Trimming preview (--trim-preview, trimming_preview.txt) on generated reads:
# 20000 reads of 100 bases whose qualities fall towards the 3' end, a third
# of them with the Illumina TruSeq adapter after an insert of 10 to 90 bases.
# The quality rows (cutoffs 10, 20, 30) equal those counted from cutadapt
# 5.2 -q on the same reads (expected_quality.tsv; with cutadapt on PATH or
# CUTADAPT=... also run directly); the same reads in Phred+64 give the same
# rows; the adapter is named; 1 and 8 threads write the same file; the JSON
# summary has the section; SAM/BAM input is refused.
# Usage: test/trim_preview/run_trim_preview_tests.sh <falcoAmadeus binary>
set -u
FA=${1:-./falcoAmadeus}
HERE=$(cd "$(dirname "$0")" && pwd)
W=$(mktemp -d)
trap 'rm -rf "$W"' EXIT
fail=0
# Park-Miller generator (awk's own rand() is not repeatable in mawk)
awk -v out="$W/r.fq" -v out64="$W/r64.fq" '
function rnd() { x = (x * 16807) % 2147483647; return x / 2147483647 }
BEGIN {
  x = 23; b = "ACGT"; L = 100
  a = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTGCTTG"
  for (k = 33; k < 127; k++) chr[k] = sprintf("%c", k)
  for (i = 0; i < 20000; i++) {
    s = ""
    for (k = 0; k < L; k++) s = s substr(b, int(rnd() * 4) + 1, 1)
    if (i % 3 == 1)
      s = substr(substr(s, 1, 10 + int(rnd() * 81)) a s, 1, L)
    q = ""; q64 = ""
    drop = 40 + int(rnd() * 60)
    for (k = 0; k < L; k++) {
      p = (k < drop) ? 30 + int(rnd() * 12) : 2 + int(rnd() * 30)
      q = q chr[p + 33]; q64 = q64 chr[p + 64]
    }
    printf "@r%d\n%s\n+\n%s\n", i, s, q > out
    printf "@r%d\n%s\n+\n%s\n", i, s, q64 > out64
  }
}'
rows() { awk -F'\t' '/^#Cutoff/ {t = 1; next} t && NF == 6 {print} t && NF != 6 {t = 0}' "$1"; }
"$FA" -q -t 1 --trim-preview -o "$W/o1" "$W/r.fq" > /dev/null 2>&1
"$FA" -q -t 8 --trim-preview -o "$W/o8" "$W/r.fq" > /dev/null 2>&1
"$FA" -q --trim-preview -o "$W/o64" "$W/r64.fq" > /dev/null 2>&1
R="$W/o1/trimming_preview.txt"
rows "$R" > "$W/ours.tsv"
if [ -s "$W/ours.tsv" ] && cmp -s "$W/ours.tsv" "$HERE/expected_quality.tsv"; then
  echo "ok   quality rows equal cutadapt -q (stored)"
else echo "FAIL quality rows differ from the stored cutadapt rows"; fail=1; fi
CA=${CUTADAPT:-$(command -v cutadapt || true)}
if [ -n "$CA" ]; then
  for c in 10 20 30; do
    "$CA" -q $c -o "$W/t$c.fq" "$W/r.fq" > /dev/null 2>&1
    python3 -I -c '
import sys
def seqs(p):
    with open(p) as f:
        return [l.rstrip("\n") for i, l in enumerate(f) if i % 4 == 1]
a, b = seqs(sys.argv[1]), seqs(sys.argv[2])
n, bases = len(a), sum(map(len, a))
rem = sum(len(x) - len(y) for x, y in zip(a, b))
tr = sum(len(y) < len(x) for x, y in zip(a, b))
sh = [sum(len(y) < m for y in b) for m in (20, 36, 50)]
g = lambda v: "%g" % v
print("\t".join([sys.argv[3], g(100.0 * rem / bases), g(100.0 * tr / n)] +
                [g(100.0 * s / n) for s in sh]))' "$W/r.fq" "$W/t$c.fq" $c
  done > "$W/ca.tsv"
  if cmp -s "$W/ours.tsv" "$W/ca.tsv"; then echo "ok   quality rows equal cutadapt -q (run)"
  else echo "FAIL quality rows differ from cutadapt"; fail=1; fi
  # REGEN_QUAL=<path>: keep cutadapt's rows (to regenerate expected_quality.tsv)
  [ -n "${REGEN_QUAL:-}" ] && cp "$W/ca.tsv" "$REGEN_QUAL"
fi
rows "$W/o64/trimming_preview.txt" > "$W/ours64.tsv"
cmp -s "$W/ours.tsv" "$W/ours64.tsv" && echo "ok   Phred+64 reads give the same rows" ||
  { echo "FAIL Phred+64 rows differ"; fail=1; }
grep -q "^Adapter	Illumina Universal Adapter" "$R" && echo "ok   adapter named" ||
  { echo "FAIL adapter: $(grep '^Adapter' "$R")"; fail=1; }
cmp -s "$R" "$W/o8/trimming_preview.txt" && echo "ok   1 and 8 threads equal" ||
  { echo "FAIL 1 and 8 threads differ"; fail=1; }
python3 -I -c 'import json,sys
d = json.load(open(sys.argv[1]))["trimming_preview"]["reads"]
assert d["quality"]["20"]["bases_removed_pct"] >= 0 and d["adapter"]' \
  "$W/o1/falcoAmadeus_summary.json" 2> /dev/null &&
  echo "ok   JSON section" || { echo "FAIL JSON section"; fail=1; }
printf '@HD\tVN:1.6\n' > "$W/x.sam"
if "$FA" -q --trim-preview -o "$W/osam" "$W/x.sam" > /dev/null 2>&1; then
  echo "FAIL SAM input accepted"; fail=1
else echo "ok   SAM input refused"; fi
exit $fail
