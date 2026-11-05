#!/usr/bin/env bash
# Per base sequence content under the methylation profiles on generated
# reads of converted composition (A 31 %, G 19 %, T 48 %, C 2 % at every
# position): the profile compares each position with the library's own mean
# composition, so an even library passes and a G block at positions 21-30 in
# half of the reads fails (RRBS composition is not graded: it passes). An
# unconverted library (all bases 25 %) is flagged by the conversion estimate.
# Usage: test/assay/run_methylation_tests.sh <falcoAmadeus binary>
set -u
FA=${1:-./falcoAmadeus}
W=$(mktemp -d)
trap 'rm -rf "$W"' EXIT
fail=0

# gen <name> <percent of reads with the G block>
# 20000 reads of 76 bases; at each position the reads cycle through a
# 100-letter pool, so every position has exactly the pool's composition
gen() {
  awk -v pct="$2" -v unconv="${3:-0}" 'BEGIN {
    pool = ""
    na = unconv ? 25 : 31; ng = unconv ? 25 : 19
    nt = unconv ? 25 : 48; nc = unconv ? 25 : 2
    for (k = 0; k < na; k++) pool = pool "A"
    for (k = 0; k < ng; k++) pool = pool "G"
    for (k = 0; k < nt; k++) pool = pool "T"
    for (k = 0; k < nc; k++) pool = pool "C"
    q = "#"; for (p = 2; p <= 76; p++) q = q "I"
    for (i = 0; i < 20000; i++) {
      s = ""
      for (p = 1; p <= 76; p++) {
        if (i % 100 < pct && p >= 21 && p <= 30) s = s "G"
        else s = s substr(pool, (i * 37 + p * 61) % 100 + 1, 1)
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
  got=$(awk -F'\t' '$2 == "Per base sequence content" {print $1}' \
          "$W/out/summary.txt")
  if [ "$got" = "$want" ]; then echo "ok   $name $*: $want"
  else echo "FAIL $name $*: '$got', expected $want"; fail=1; fi
}

gen even 0
gen gblock 50
gen unconv 0 1

for a in wgbs emseq; do
  expect even PASS --assay "$a"
  expect gblock FAIL --assay "$a"
done
expect even PASS --assay rrbs
expect gblock PASS --assay rrbs  # RRBS composition is not graded
expect even WARN  # generic grading: |A-T| and |C-G| are 17
# unconverted reads: even composition, but the conversion estimate (C among
# C+T outside CpG, about 50 %) is above the 20 % limit
expect unconv WARN --assay wgbs
if grep -q "possible incomplete bisulfite conversion" "$W/out/fastqc_data.txt"
then echo "ok   unconv --assay wgbs: conversion warning"
else echo "FAIL unconv --assay wgbs: no conversion warning"; fail=1; fi

exit $fail
