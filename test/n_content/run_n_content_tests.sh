#!/usr/bin/env bash
# Per base N content grades on generated reads with known N fractions.
# By default the module is graded as in FastQC (warn above 5 %, fail above
# 20 % N per position or position group); --legacy-n-content-grade keeps the
# falco 1.x grade, which never warns or fails.
# Usage: test/n_content/run_n_content_tests.sh <falcoAmadeus binary>
set -u
FA=${1:-./falcoAmadeus}
W=$(mktemp -d)
trap 'rm -rf "$W"' EXIT
fail=0

# gen <name> <percent of reads with N> <first N position> <last N position>
# 20000 reads of 76 bases; reads i with i % 100 < percent carry N at the
# given 1-based positions
gen() {
  awk -v pct="$2" -v a="$3" -v b="$4" 'BEGIN {
    q = "#"; for (p = 2; p <= 76; p++) q = q "I"
    for (i = 0; i < 20000; i++) {
      s = ""
      for (p = 1; p <= 76; p++) {
        if (i % 100 < pct && p >= a && p <= b) s = s "N"
        else s = s substr("ACGT", (i * 7 + p * 13) % 4 + 1, 1)
      }
      printf "@r%d\n%s\n+\n%s\n", i, s, q
    }
  }' > "$W/$1.fastq"
}

# expect <name> <grade> [options...]
expect() {
  local name=$1 want=$2 got
  shift 2
  rm -rf "$W/out"
  mkdir -p "$W/out"
  "$FA" -q "$@" -o "$W/out" "$W/$name.fastq" > /dev/null 2>&1
  got=$(awk -F'\t' '$2 == "Per base N content" {print $1}' "$W/out/summary.txt")
  if [ "$got" = "$want" ]; then echo "ok   $name $*: $want"
  else echo "FAIL $name $*: '$got', expected $want"; fail=1; fi
}

gen n03 3 11 20
gen n10 10 11 20
gen n50 50 11 20
gen n1pos 30 40 40  # one position: 6 % over its group of five bases

expect n03 PASS
expect n10 WARN
expect n50 FAIL
expect n1pos WARN
expect n1pos FAIL --nogroup
for n in n03 n10 n50 n1pos; do
  expect $n PASS --legacy-n-content-grade
done

exit $fail
