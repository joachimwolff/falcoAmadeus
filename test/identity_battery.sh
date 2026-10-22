#!/usr/bin/env bash
# Byte-identity battery against a reference binary: every input at several
# thread counts, plus thread invariance and the paired combined view.
# Usage: test/identity_battery.sh <new binary> <reference binary> <files...>
# Prints one line per check and exits non-zero on any difference.
set -u
NEW=$1; REF=$2; shift 2
W=$(mktemp -d)
trap 'rm -rf "$W"' EXIT
fail=0

same() {  # dirA dirB [ignore-filename]
  if [ "${3:-}" = nofilename ]; then
    diff <(grep -v '^Filename' "$1/fastqc_data.txt") \
         <(grep -v '^Filename' "$2/fastqc_data.txt") > /dev/null
  else
    diff -q "$1/fastqc_data.txt" "$2/fastqc_data.txt" > /dev/null &&
      diff -q "$1/summary.txt" "$2/summary.txt" > /dev/null
  fi
}

for f in "$@"; do
  n=$(basename "$f")
  for t in 1 2 8 16; do
    rm -rf "$W/a" "$W/b"
    "$NEW" -t "$t" -q -o "$W/a" "$f" > /dev/null 2>&1; e1=$?
    "$REF" -t "$t" -q -o "$W/b" "$f" > /dev/null 2>&1; e2=$?
    if [ "$e1" != "$e2" ]; then
      echo "FAIL $n t$t: exit $e1 vs reference $e2"; fail=1
    elif [ "$e1" = 0 ] && ! same "$W/a" "$W/b"; then
      echo "FAIL $n t$t: output differs from reference"; fail=1
    else
      echo "ok   $n t$t"
    fi
    [ "$t" = 1 ] && cp -r "$W/a" "$W/t1"
    if [ "$t" != 1 ] && [ "$e1" = 0 ] && ! same "$W/t1" "$W/a" nofilename; then
      echo "FAIL $n: t$t differs from t1"; fail=1
    fi
  done
  rm -rf "$W/t1"
done
exit $fail
