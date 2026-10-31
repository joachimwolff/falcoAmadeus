#!/usr/bin/env bash
# Library complexity (--complexity) against preseq 3.2.0 lc_extrap -Q: a
# generated library (20000 molecules, copy numbers drawn deterministically
# from a skewed distribution) gives the same yield curve as preseq on the
# same copy-number histogram. expected_curve.tsv holds preseq's output; with
# preseq on PATH (or PRESEQ=...) it is also run directly.
# Usage: test/complexity/run_complexity_tests.sh <falcoAmadeus binary>
set -u
FA=${1:-./falcoAmadeus}
HERE=$(cd "$(dirname "$0")" && pwd)
W=$(mktemp -d)
trap 'rm -rf "$W"' EXIT
fail=0
# a Park-Miller generator: awk's srand/rand differ between implementations
# (mawk does not repeat a seeded sequence); the products stay below 2^46,
# exact in double precision
awk 'function rnd() { x = (x * 16807) % 2147483647; return x / 2147483647 }
BEGIN {
  x = 7; b = "ACGT"
  q = "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII"
  for (m = 0; m < 20000; m++) {
    s = ""
    for (k = 0; k < 60; k++) s = s substr(b, int(rnd() * 4) + 1, 1)
    c = 1 + int(-log(1 - rnd()) * (m % 7 == 0 ? 6 : 1.2))
    for (j = 0; j < c; j++) printf "@m%d_%d\n%s\n+\n%s\n", m, j, s, q
  }
}' > "$W/lib.fastq"
FA_COPY_HIST="$W/hist.txt" FA_COMPLEXITY_CURVE="$W/curve.txt" \
  "$FA" -q --complexity -o "$W/out" "$W/lib.fastq" > /dev/null 2>&1
if cmp -s "$W/curve.txt" "$HERE/expected_curve.tsv"; then
  echo "ok   curve equals preseq lc_extrap -Q (stored)"
else
  echo "FAIL curve differs from the stored preseq output"; fail=1
fi
PQ=${PRESEQ:-$(command -v preseq || true)}
if [ -n "$PQ" ]; then
  "$PQ" lc_extrap -Q -H "$W/hist.txt" -o "$W/pq.txt" 2> /dev/null
  if cmp -s "$W/curve.txt" "$W/pq.txt"; then echo "ok   curve equals preseq lc_extrap -Q (run)"
  else echo "FAIL curve differs from preseq"; fail=1; fi
fi
grep -q "^10x" "$W/out/complexity.txt" && echo "ok   complexity.txt written" ||
  { echo "FAIL complexity.txt"; fail=1; }
# REGEN_HIST=<path>: keep the histogram (to regenerate expected_curve.tsv)
[ -n "${REGEN_HIST:-}" ] && cp "$W/hist.txt" "$REGEN_HIST"
exit $fail
