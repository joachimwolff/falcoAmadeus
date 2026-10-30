#!/usr/bin/env bash
# Error- and UMI-aware paired duplication against an all-pairs oracle
# (aware_oracle.py) on simulated read pairs, paired and single-end.
# Usage: run_aware_tests.sh <falcoAmadeus>
set -euo pipefail
FA=$1
HERE=$(cd "$(dirname "$0")" && pwd)
T=$(mktemp -d)
trap 'rm -rf "$T"' EXIT
fail=0
for seed in 1 2; do
  python3 -I "$HERE/aware_oracle.py" simulate $seed "$T/r1.fq" "$T/r2.fq"
  while read -r D U UMI; do
    args=(--dup-mismatches "$D")
    case $UMI in
      name) args+=(--umi-in-name) ;;
      r1:*) args+=(--umi-r1 "${UMI#r1:}") ;;
      r1r2:*) args+=(--umi-r1 "${UMI#r1r2:}" --umi-r2 "${UMI#r1r2:}") ;;
    esac
    [ "$U" != 0 ] && args+=(--umi-mismatches "$U")
    for t in 1 4; do
      rm -rf "$T/out"
      FA_AWARE_CLUSTERS="$T/g.u32" "$FA" -q -t $t -o "$T/out" "$T/r1.fq" \
        --paired-with "$T/r2.fq" "${args[@]}" > /dev/null
      python3 -I "$HERE/aware_oracle.py" check "$T/r1.fq" "$T/r2.fq" \
        "$T/g.u32" "$D" "$U" "$UMI" || fail=1
      # single-end: read 1 alone (no UMI in read 2)
      case $UMI in r1r2:*) continue ;; esac
      rm -rf "$T/out"
      FA_AWARE_CLUSTERS="$T/g.u32" "$FA" -q -t $t -o "$T/out" "$T/r1.fq" \
        "${args[@]}" > /dev/null
      [ -s "$T/out/aware_duplication.txt" ] || { echo "FAIL no SE report"; fail=1; }
      python3 -I "$HERE/aware_oracle.py" check "$T/r1.fq" - \
        "$T/g.u32" "$D" "$U" "$UMI" || fail=1
    done
  done <<'CASES'
1 0 none
2 0 none
3 0 none
5 0 none
0 0 name
0 1 name
2 1 name
0 0 r1:1:6
1 1 r1:1:6
2 0 r1r2:1:6
2 2 r1r2:1:6
CASES
done
exit $fail
