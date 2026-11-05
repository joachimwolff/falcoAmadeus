#!/usr/bin/env bash
# Paired-end duplicate rates on simulated read pairs with a known duplicate
# fraction (test/paired/sim_paired_dup.py): the exact and FastQC-style rates
# reported by falcoAmadeus --paired-with, for three placements of the copies,
# and for copies that repeat R1 only (true pair duplicate rate 0).
# Usage: run_duplication_simulation.sh <falcoAmadeus> <out.tsv> [pairs]
set -eu
FA=$1; OUT=$2; PAIRS=${3:-300000}
HERE=$(cd "$(dirname "$0")" && pwd)
SIM=$HERE/../../test/paired/sim_paired_dup.py
W=$(mktemp -d)
trap 'rm -rf "$W"' EXIT
get() { awk -F'\t' -v k="$2" '$1 == k {print $2}' "$1"; }
echo -e "scenario\tplacement\ttrue_dup\tpaired_exact\tr1_exact\tr2_exact\tpaired_fastqc\tr1_fastqc\tr2_fastqc" > "$OUT"
run() {  # scenario placement dup [sim flags]
  local sc=$1 pl=$2 d=$3; shift 3
  python3 -I "$SIM" --prefix "$W/s" --pairs "$PAIRS" --dup "$d" --placement "$pl" --seed 7 "$@"
  rm -rf "$W/o"; mkdir "$W/o"
  "$FA" -t 4 -q -o "$W/o" "$W/s_R1.fastq" --paired-with "$W/s_R2.fastq" > /dev/null
  local f="$W/o/paired_duplication.txt"
  echo -e "$sc\t$pl\t$d\t$(get $f 'Paired-end exact')\t$(get $f 'Read1-only exact')\t$(get $f 'Read2-only exact')\t$(get $f 'Paired-end FastQC-style')\t$(get $f 'Read1-only FastQC-style')\t$(get $f 'Read2-only FastQC-style')" >> "$OUT"
}
for d in 0 0.1 0.2 0.3 0.45 0.6; do
  for pl in late early shuffled; do run both_mates $pl $d; done
  run shared_r1 shuffled $d --shared-r1
done
