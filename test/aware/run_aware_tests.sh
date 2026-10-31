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
# a UMI that ends the name line (no comment after it, a /1 /2 suffix)
python3 -I - "$T" <<'PY'
import sys
d = sys.argv[1]
with open(d + '/e1.fq', 'w') as a, open(d + '/e2.fq', 'w') as b:
    for k in range(20):
        u = 'ACGTACGT' if k % 2 else 'TTGCAAGC'
        a.write('@r%d:%s/1\n%s\n+\n%s\n' % (k, u, 'ACGT' * 25, 'I' * 100))
        b.write('@r%d:%s/2\n%s\n+\n%s\n' % (k, u, 'TGCA' * 25, 'I' * 100))
PY
rm -rf "$T/out"
if "$FA" -q -o "$T/out" "$T/e1.fq" --paired-with "$T/e2.fq" --umi-in-name \
    > /dev/null 2>&1 &&
  grep -q "^Pairs without a UMI in the read name	0$" "$T/out/paired_duplication.txt" &&
  grep -q "^Groups of copies	2$" "$T/out/paired_duplication.txt"; then
  echo "ok   UMI at the end of the name line"
else echo "FAIL UMI at the end of the name line"; fail=1; fi
exit $fail
