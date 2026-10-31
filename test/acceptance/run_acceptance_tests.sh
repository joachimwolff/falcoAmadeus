#!/usr/bin/env bash
# Acceptance rules (--acceptance, --acceptance-exit) on a small library:
# each rule graded PASS/WARN/FAIL/MISSING, the overall result, exit status 2
# on FAIL with --acceptance-exit (0 without it), and a malformed rules file
# rejected before the analysis.
# Usage: test/acceptance/run_acceptance_tests.sh <falcoAmadeus binary>
set -u
FA=${1:-./falcoAmadeus}
W=$(mktemp -d)
trap 'rm -rf "$W"' EXIT
fail=0
awk 'BEGIN {
  x = 9
  for (i = 0; i < 2000; i++) {
    s = ""; q = "#"  # one Phred-2 base: Sanger encoding (else Illumina 1.5)
    for (k = 0; k < 60; k++) {
      x = (x * 16807) % 2147483647
      s = s substr("ACGT", int(x / 2147483647 * 4) + 1, 1)
      if (k > 0) q = q "I"
    }
    printf "@r%d\n%s\n+\n%s\n", i, s, q
  }
}' > "$W/r.fq"
cat > "$W/rules" <<'R'
q30_bases_pct                                >= 98  90    # 59 of 60 bases Q40
views.reads.basic_statistics.total_sequences >= 5000 1000 # 2000 reads: WARN
views.reads.basic_statistics.total_sequences <= 1500 1800 # above 1800: FAIL
no.such.metric                               >= 1   0
R
"$FA" -q --acceptance "$W/rules" -o "$W/o" "$W/r.fq" > /dev/null 2>&1
st=$?
A="$W/o/acceptance.txt"
expect() {  # line-number result
  got=$(awk -F'\t' -v n="$1" '/^#Metric/ {f = 1; k = 0; next} f {k++; if (k == n) print $6}' "$A")
  [ "$got" = "$2" ] && echo "ok   rule $1: $2" || { echo "FAIL rule $1: $got, expected $2"; fail=1; }
}
expect 1 PASS
expect 2 WARN
expect 3 FAIL
expect 4 MISSING
grep -q "^Overall	FAIL" "$A" && echo "ok   overall FAIL" || { echo "FAIL overall"; fail=1; }
[ $st = 0 ] && echo "ok   exit 0 without --acceptance-exit" || { echo "FAIL exit $st"; fail=1; }
"$FA" -q --acceptance "$W/rules" --acceptance-exit -o "$W/o2" "$W/r.fq" > /dev/null 2>&1
st=$?
[ $st = 2 ] && echo "ok   exit 2 with --acceptance-exit" || { echo "FAIL exit $st, expected 2"; fail=1; }
printf 'q30_bases_pct >= ninety 80\n' > "$W/bad"
"$FA" -q --acceptance "$W/bad" -o "$W/o3" "$W/r.fq" > /dev/null 2>&1
st=$?
[ $st = 1 ] && [ ! -e "$W/o3/fastqc_data.txt" ] &&
  echo "ok   malformed rules rejected before the analysis" ||
  { echo "FAIL malformed rules: exit $st"; fail=1; }
exit $fail
