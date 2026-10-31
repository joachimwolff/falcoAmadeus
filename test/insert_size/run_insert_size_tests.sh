#!/usr/bin/env bash
# Insert size from the mate overlap (paired_insert_size.txt) on generated
# pairs: fragments of 40 to 400 bases, 100-base reads, the Illumina TruSeq
# adapters after fragments shorter than a read, and 400 adapter dimers in
# 20000 pairs. The histogram equals fastp 0.24.3's on the same pairs
# (expected_hist.tsv; with fastp on PATH or FASTP=... also run directly),
# the adapter of read 1 is named, and every dimer is found.
# Usage: test/insert_size/run_insert_size_tests.sh <falcoAmadeus binary>
set -u
FA=${1:-./falcoAmadeus}
HERE=$(cd "$(dirname "$0")" && pwd)
W=$(mktemp -d)
trap 'rm -rf "$W"' EXIT
fail=0
# Park-Miller generator (awk's own rand() is not repeatable in mawk)
awk -v out1="$W/r1.fq" -v out2="$W/r2.fq" '
function rnd() { x = (x * 16807) % 2147483647; return x / 2147483647 }
function rc(s,   r, i, c) {
  r = ""
  for (i = length(s); i >= 1; i--) {
    c = substr(s, i, 1)
    r = r (c == "A" ? "T" : c == "C" ? "G" : c == "G" ? "C" : "A")
  }
  return r
}
BEGIN {
  x = 11; b = "ACGT"; L = 100
  a1 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTGCTTGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
  a2 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATTAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
  q = ""; for (k = 0; k < L; k++) q = q "I"
  for (i = 0; i < 20000; i++) {
    if (i % 50 == 7) f = ""  # adapter dimer: no insert
    else {
      n = 40 + int(rnd() * 361); f = ""
      for (k = 0; k < n; k++) f = f substr(b, int(rnd() * 4) + 1, 1)
    }
    s1 = substr(f a1, 1, L); s2 = substr(rc(f) a2, 1, L)
    printf "@p%d/1\n%s\n+\n%s\n", i, s1, q > out1
    printf "@p%d/2\n%s\n+\n%s\n", i, s2, q > out2
  }
}'
"$FA" -q -o "$W/out" "$W/r1.fq" --paired-with "$W/r2.fq" > /dev/null 2>&1
R="$W/out/paired_insert_size.txt"
awk -F'\t' '$1 ~ /^[0-9]+$/ && NF == 2 {h[$1] = $2} $1 ~ /^>=512/ {h[512] = $2}
  END {for (i = 0; i <= 512; i++) print (i in h) ? h[i] : 0}' "$R" > "$W/ours.txt"
if cmp -s "$W/ours.txt" "$HERE/expected_hist.tsv"; then
  echo "ok   histogram equals fastp (stored)"
else echo "FAIL histogram differs from the stored fastp histogram"; fail=1; fi
FP=${FASTP:-$(command -v fastp || true)}
if [ -n "$FP" ]; then
  "$FP" -w 1 -A -G -i "$W/r1.fq" -I "$W/r2.fq" -o /dev/null -O "$W/o2.fq" \
    -j "$W/fp.json" -h /dev/null 2> /dev/null
  python3 -I -c 'import json,sys; j=json.load(open(sys.argv[1]))["insert_size"]
for v in j["histogram"] + [j["unknown"]]: print(v)' "$W/fp.json" > "$W/fp.txt"
  if cmp -s "$W/ours.txt" "$W/fp.txt"; then echo "ok   histogram equals fastp (run)"
  else echo "FAIL histogram differs from fastp"; fail=1; fi
fi
grep -q "^Adapter read 1.*Illumina Universal Adapter" "$R" &&
  echo "ok   adapter of read 1 named" ||
  { echo "FAIL adapter of read 1: $(grep '^Adapter read 1' "$R")"; fail=1; }
d=$(awk -F'\t' '$1 == "Adapter dimers" {print $2}' "$R")
[ "$d" = 400 ] && echo "ok   400 adapter dimers found" ||
  { echo "FAIL adapter dimers: $d, expected 400"; fail=1; }
# REGEN_HIST=<path>: keep our histogram (to regenerate expected_hist.tsv)
[ -n "${REGEN_HIST:-}" ] && cp "$W/ours.txt" "$REGEN_HIST"
exit $fail
