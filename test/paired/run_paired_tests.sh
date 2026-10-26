#!/usr/bin/env bash
# Paired-end duplication and pairing tests on simulated data with known truth.
# Usage: test/paired/run_paired_tests.sh <falcoAmadeus binary>
set -u
FA=${1:-./falcoAmadeus}
HERE=$(cd "$(dirname "$0")" && pwd)
W=$(mktemp -d)
trap 'rm -rf "$W"' EXIT
fail=0

sim() { python3 -I "$HERE/sim_paired_dup.py" --prefix "$W/$1" "${@:2}"; }

# expect <name> <key> <value>: exact match of a paired_duplication.txt row
expect() {
  local got
  got=$(awk -F'\t' -v k="$2" '$1 == k {print $2}' "$W/out_$1/paired_duplication.txt")
  if [ "$got" = "$3" ]; then echo "ok   $1: $2 = $3"
  else echo "FAIL $1: $2 = '$got', expected $3"; fail=1; fi
}

run() {
  mkdir -p "$W/out_$1"
  "$FA" -t 4 -q -o "$W/out_$1" "$W/${1}_R1.fastq" \
    --paired-with "$W/${1}_R2.fastq" > "$W/$1.log" 2>&1
}

sim unique --pairs 1000 --dup 0 && run unique
expect unique "Paired-end exact" 0

sim dup25 --pairs 1000 --dup 0.25 && run dup25
expect dup25 "Paired-end exact" 25

# copies read from the other strand (mates swapped) are the same fragment
sim swap25 --pairs 1000 --dup 0.25 --swap-mates && run swap25
expect swap25 "Paired-end exact" 25

# exact rate must not depend on where the copies sit in the file
for p in late early shuffled; do
  sim "o_$p" --pairs 200000 --dup 0.45 --placement "$p" && run "o_$p"
  expect "o_$p" "Paired-end exact" 45
done

# copies share R1 only: paired 0, R1-only 25
sim shared --pairs 1000 --dup 0.25 --shared-r1 && run shared
expect shared "Paired-end exact" 0
expect shared "Read1-only exact" 25

# mates out of order must be rejected
sim reorder --pairs 1000 --dup 0.25 --reorder-r2
if run reorder; then echo "FAIL reorder: accepted mismatched mates"; fail=1
else echo "ok   reorder: rejected mismatched mates"; fi

# fastq-dump --readids naming (<spot>.1 / <spot>.2) pairs; shifted mates do not
awk 'NR % 4 == 1 {$1 = $1 ".1"} {print}' "$W/dup25_R1.fastq" > "$W/sra_R1.fastq"
awk 'NR % 4 == 1 {$1 = $1 ".2"} {print}' "$W/dup25_R2.fastq" > "$W/sra_R2.fastq"
run sra
expect sra "Paired-end exact" 25
tail -n +5 "$W/sra_R2.fastq" > "$W/srashift_R2.fastq"
head -n -4 "$W/sra_R1.fastq" > "$W/srashift_R1.fastq"
if run srashift; then echo "FAIL srashift: accepted shifted mates"; fail=1
else echo "ok   srashift: rejected shifted mates"; fi

# the combined view equals the analysis of R1 followed by R2 (with more than
# 100k distinct R1 sequences it is derived from the two mates)
sim pool --pairs 150000 --dup 0.2 --placement shuffled
for t in 1 4; do
  mkdir -p "$W/pool$t" "$W/cat$t"
  "$FA" -t "$t" -q --paired-data both -o "$W/pool$t" "$W/pool_R1.fastq" \
    --paired-with "$W/pool_R2.fastq" > /dev/null 2>&1
  cat "$W/pool_R1.fastq" "$W/pool_R2.fastq" > "$W/pool_cat.fastq"
  "$FA" -t "$t" -q -o "$W/cat$t" "$W/pool_cat.fastq" > /dev/null 2>&1
  if diff <(grep -v '^Filename' "$W/pool$t/fastqc_data.txt") \
          <(grep -v '^Filename' "$W/cat$t/fastqc_data.txt") > /dev/null; then
    echo "ok   combined view = analysis of R1+R2 (-t $t)"
  else echo "FAIL combined view differs from R1+R2 (-t $t)"; fail=1; fi
done

# thread invariance of the paired report
for t in 1 8; do
  mkdir -p "$W/t$t"
  "$FA" -t "$t" -q -o "$W/t$t" "$W/dup25_R1.fastq" \
    --paired-with "$W/dup25_R2.fastq" > /dev/null 2>&1
done
if diff <(tail -n +4 "$W/t1/paired_duplication.txt") \
        <(tail -n +4 "$W/t8/paired_duplication.txt") > /dev/null; then
  echo "ok   paired report identical at -t 1 and -t 8"
else echo "FAIL paired report differs between -t 1 and -t 8"; fail=1; fi

exit $fail
