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
# short inserts followed by a shared adapter: unrelated pairs share their
# adapter blocks (block groups above the refinement size), and copies carry
# an error in every block of the insert, so the adapter blocks are the only
# blocks they share (single-end with 2 mismatches: blocks of 17 bases;
# paired with 4: blocks of 10). Found only by the refinement of large block
# groups.
python3 -I - "$T" <<'PY'
import random, sys
d = sys.argv[1]
A1 = 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATC'
A2 = 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGT'
def gen(seed, out1, out2, ins_len, cuts):
    rnd = random.Random(seed)
    def hit(s, lo, hi):
        s = list(s)
        i = rnd.randrange(lo, hi)
        s[i] = rnd.choice('ACGT'.replace(s[i], ''))
        return ''.join(s)
    f1, f2 = open(out1, 'w'), (open(out2, 'w') if out2 else None)
    k = 0
    for m in range(600):
        ins = ''.join(rnd.choice('ACGT') for _ in range(rnd.randrange(*ins_len)))
        rc = ins[::-1].translate(str.maketrans('ACGT', 'TGCA'))
        a0, b0 = (ins + A1)[:51], (rc + A2)[:51]
        for c in range(1 + (m % 3)):
            a, b = a0, b0
            if c:
                for lo, hi in cuts:
                    a = hit(a, lo, hi)
                    if out2:
                        b = hit(b, lo, hi)
            if c == 2 and out2:
                a, b = b, a
            f1.write('@p%d\n%s\n+\n%s\n' % (k, a, 'I' * 51))
            if out2:
                f2.write('@p%d\n%s\n+\n%s\n' % (k, b, 'I' * 51))
            k += 1
gen(11, d + '/se.fq', None, (18, 31), [(0, 17), (17, 34)])
gen(12, d + '/pe1.fq', d + '/pe2.fq', (14, 21), [(0, 10), (10, 20)])
PY
# FA_REFINE_MIN 2: every group above two classes split, as deep as it goes
for rmin in 64 2; do
  rm -rf "$T/out"
  FA_REFINE_MIN=$rmin FA_AWARE_CLUSTERS="$T/g.u32" "$FA" -q -t 4 -o "$T/out" \
    "$T/se.fq" --dup-mismatches 2 > /dev/null
  python3 -I "$HERE/aware_oracle.py" check "$T/se.fq" - "$T/g.u32" 2 0 none || fail=1
  rm -rf "$T/out"
  FA_REFINE_MIN=$rmin FA_AWARE_CLUSTERS="$T/g.u32" "$FA" -q -t 4 -o "$T/out" \
    "$T/pe1.fq" --paired-with "$T/pe2.fq" --dup-mismatches 4 > /dev/null
  python3 -I "$HERE/aware_oracle.py" check "$T/pe1.fq" "$T/pe2.fq" "$T/g.u32" 4 0 none || fail=1
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
