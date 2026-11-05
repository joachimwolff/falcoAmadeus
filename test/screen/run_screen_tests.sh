#!/usr/bin/env bash
# Contamination screen (--screen-build, --screen) on generated references:
# genome A (400 kb), B (80 kb), C (5 kb, kept whole) and reads mixed 80 % A,
# 15 % B, 4 % C and 1 % low-complexity (poly-A and CA repeats). Without
# sequencing errors each share of found k-mers must be within 10 % (relative)
# of the read share, the low-complexity reads must match nothing, and a
# database must give the same result after saving and loading.
# Usage: test/screen/run_screen_tests.sh <falcoAmadeus binary>
set -u
FA=${1:-./falcoAmadeus}
W=$(mktemp -d)
trap 'rm -rf "$W"' EXIT
fail=0
# Park-Miller generator (awk's own rand() is not repeatable in mawk)
awk -v d="$W" '
function rnd() { x = (x * 16807) % 2147483647; return x / 2147483647 }
function genome(name, n,   s, k, j, c) {
  s = ""
  for (k = 0; k < n; k += 1000) {  # in chunks: appending base by base is
    c = ""                          # quadratic
    for (j = 0; j < 1000; j++) c = c substr("ACGT", int(rnd() * 4) + 1, 1)
    s = s c
  }
  printf ">%s\n%s\n", name, s > (d "/" name ".fa")
  close(d "/" name ".fa")
  return s
}
BEGIN {
  x = 3
  A = genome("A", 400000); B = genome("B", 80000); C = genome("C", 5000)
  q = ""; for (k = 0; k < 100; k++) q = q "I"
  lc = ""; for (k = 0; k < 50; k++) lc = lc "A"
  for (k = 0; k < 25; k++) lc = lc "CA"
  for (i = 0; i < 100000; i++) {
    r = rnd()
    if (r < 0.80) g = A; else if (r < 0.95) g = B; else if (r < 0.99) g = C
    else g = ""
    s = (g == "") ? lc : substr(g, int(rnd() * (length(g) - 100)) + 1, 100)
    printf "@r%d\n%s\n+\n%s\n", i, s, q > (d "/reads.fq")
  }
}'
"$FA" -q -t 2 --screen-build "$W/db.fascreen" \
  --screen-refs "A=$W/A.fa,B=$W/B.fa,C=$W/C.fa" > /dev/null 2>&1 ||
  { echo "FAIL --screen-build"; exit 1; }
"$FA" -q -t 2 --screen "$W/db.fascreen" -o "$W/out" "$W/reads.fq" \
  > /dev/null 2>&1 || { echo "FAIL --screen"; exit 1; }
R="$W/out/screen.txt"
check() {  # name expected-percent
  awk -F'\t' -v n="$1" -v e="$2" '$1 == n {
    d = ($2 - e) / e; if (d < 0) d = -d
    if (d <= 0.10) printf "ok   %s found %.2f %% (expected about %s)\n", n, $2, e
    else { printf "FAIL %s found %.2f %%, expected about %s\n", n, $2, e; bad = 1 }
  } END {exit bad}' "$R" || fail=1
}
# expected: read share of each genome's k-mers among all k-mers that pass the
# low-complexity filter (the low-complexity reads contribute none)
check A 80.8
check B 15.2
check C 4.0
awk -F'\t' '$1 == "No reference" && $2 < 0.5 {ok = 1} END {exit !ok}' "$R" &&
  echo "ok   low-complexity reads add no unmatched k-mers" ||
  { echo "FAIL unmatched k-mers: $(grep '^No reference' "$R")"; fail=1; }
cp "$W/db.fascreen" "$W/db2.fascreen"
"$FA" -q -t 1 --screen "$W/db2.fascreen" -o "$W/out2" "$W/reads.fq" \
  > /dev/null 2>&1
if diff <(sed 1,2d "$R") <(sed 1,2d "$W/out2/screen.txt") > /dev/null; then
  echo "ok   same result at -t 1 from a copied database"
else echo "FAIL result differs at -t 1 / copied database"; fail=1; fi
exit $fail
