#!/usr/bin/env bash
# UMI quality (umi_qc.txt) on generated reads: molecules with family sizes
# drawn from an exponential, base errors in UMIs and inserts, N and poly-G
# UMIs. Every line of umi_qc.txt (UMI sequences, families, UMI base quality)
# equals umi_oracle.py, an independent computation from the FASTQ files and
# the groups of the same run (FA_AWARE_CLUSTERS). Cases: paired with UMIs in
# both mates; single-end with reads shorter than the UMI region; UMIs of
# varying length in the read names; Phred+64 qualities; gzip input. 1 and 16
# threads write the same file; the JSON summary has the section.
# Usage: test/umi_qc/run_umi_qc_tests.sh <falcoAmadeus binary>
set -u
FA=${1:-./falcoAmadeus}
HERE=$(cd "$(dirname "$0")" && pwd)
W=$(mktemp -d)
trap 'rm -rf "$W"' EXIT
fail=0
python3 -I - "$W" <<'PY'
import random, sys, gzip
d = sys.argv[1]
B = 'ACGT'
def gen(seed, n_mol, names=False, short=False, off=33):
    rnd = random.Random(seed)
    rs = lambda n: ''.join(rnd.choice(B) for _ in range(n))
    mut = lambda s, p: ''.join(c if rnd.random() > p else
                               rnd.choice(B.replace(c, '')) for c in s)
    r1, r2 = [], []
    k = 0
    for m in range(n_mol):
        u1 = rs(8) if m % 50 else 'G' * 8
        u2 = rs(6)
        un = rs(rnd.choice((8, 10, 12)))
        i1, i2 = rs(92), rs(94)
        for c in range(1 + int(rnd.expovariate(0.4))):
            a = mut(u1, 0.01) + mut(i1, 0.003)
            b = mut(u2, 0.01) + mut(i2, 0.003)
            if rnd.random() < 0.01:
                a = 'N' + a[1:]
            if short and rnd.random() < 0.02:
                a = a[:rnd.choice((3, 6, 0))] or 'A'
            q1 = ''.join(chr(off + rnd.choice((2, 12, 25, 37, 37))) for _ in a)
            q2 = ''.join(chr(off + 37) for _ in b)
            h = '@r%d' % k + (':' + mut(un, 0.01) if names else '')
            r1.append('%s/1\n%s\n+\n%s\n' % (h, a, q1))
            r2.append('%s/2\n%s\n+\n%s\n' % (h, b, q2))
            k += 1
    return r1, r2
r1, r2 = gen(5, 3000)
open(d + '/p1.fq', 'w').write(''.join(r1))
open(d + '/p2.fq', 'w').write(''.join(r2))
s1, _ = gen(6, 3000, short=True)
open(d + '/s.fq', 'w').write(''.join(s1))
n1, n2 = gen(7, 3000, names=True)
open(d + '/n1.fq', 'w').write(''.join(n1))
open(d + '/n2.fq', 'w').write(''.join(n2))
h1, _ = gen(8, 2000, off=64)
open(d + '/h.fq', 'w').write(''.join(h1))
with gzip.open(d + '/p1.fq.gz', 'wt') as f:
    f.write(''.join(r1))
with gzip.open(d + '/p2.fq.gz', 'wt') as f:
    f.write(''.join(r2))
PY
check() {  # name, oracle inputs/options, falcoAmadeus inputs/options
  local name=$1 o=$2 f=$3
  rm -rf "$W/o" "$W/cl.bin"
  FA_AWARE_CLUSTERS="$W/cl.bin" "$FA" -q -o "$W/o" $f > /dev/null 2>&1
  local r
  r=$(ls "$W"/o/umi_qc.txt "$W"/o/*umi_qc.txt 2> /dev/null | head -1)
  if [ -z "$r" ]; then echo "FAIL $name: no umi_qc.txt"; fail=1; return; fi
  python3 -I "$HERE/umi_oracle.py" "$W/cl.bin" $o > "$W/oracle.txt"
  if diff <(sed 1,2d "$r") "$W/oracle.txt" > /dev/null; then
    echo "ok   $name: equals the oracle"
  else echo "FAIL $name: differs from the oracle"; fail=1; fi
  cp "$r" "$W/last_umi_qc.txt"
}
U="--umi-r1 1:8 --umi-r2 1:6"
check "paired, UMIs in both mates" "$W/p1.fq $W/p2.fq $U" \
  "$W/p1.fq --paired-with $W/p2.fq $U --dup-mismatches 2 --umi-mismatches 1"
check "paired, exact" "$W/p1.fq $W/p2.fq $U" "$W/p1.fq --paired-with $W/p2.fq $U"
check "single-end, short reads" "$W/s.fq --umi-r1 1:8" \
  "$W/s.fq --umi-r1 1:8 --dup-mismatches 1"
check "paired, UMIs in names" "$W/n1.fq $W/n2.fq --umi-in-name" \
  "$W/n1.fq --paired-with $W/n2.fq --umi-in-name --umi-mismatches 1"
check "single-end, Phred+64" "$W/h.fq --umi-r1 1:8" "$W/h.fq --umi-r1 1:8"
check "paired, gzip" "$W/p1.fq $W/p2.fq $U" \
  "$W/p1.fq.gz --paired-with $W/p2.fq.gz $U --dup-mismatches 2"
rm -rf "$W/t1" "$W/t16"
"$FA" -q -t 1 -o "$W/t1" "$W/p1.fq" --paired-with "$W/p2.fq" $U --dup-mismatches 2 > /dev/null 2>&1
"$FA" -q -t 16 -o "$W/t16" "$W/p1.fq" --paired-with "$W/p2.fq" $U --dup-mismatches 2 > /dev/null 2>&1
cmp -s "$W/t1/umi_qc.txt" "$W/t16/umi_qc.txt" && echo "ok   1 and 16 threads equal" ||
  { echo "FAIL 1 and 16 threads differ"; fail=1; }
python3 -I -c 'import json,sys
d = json.load(open(sys.argv[1]))["umi_qc"]
assert d["families"] > 0 and d["read_1_umi_distinct"] > 0' \
  "$W/t16/falcoAmadeus_summary.json" 2> /dev/null &&
  echo "ok   JSON section" || { echo "FAIL JSON section"; fail=1; }
exit $fail
