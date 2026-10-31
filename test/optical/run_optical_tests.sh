#!/usr/bin/env bash
# Optical duplicates (--optical) on read pairs from make_optical_data.py:
# copies are sequence-identical and at the same coordinates, so the groups
# of copies equal Picard's duplicate sets. Duplicates and optical duplicates
# equal those of Picard MarkDuplicates (GATK 4.6.2.0) on the matching
# alignments at 50, 100 and 2500 pixels (expected_picard.tsv; with GATK=gatk
# and samtools on PATH also run directly), for 1 and 16 threads. Names
# without coordinates give no optical duplicates, and single-end input is
# reported per read. Pairs that would not align (adapter dimers, inserts
# under 30 bases, mostly N) are left out, as Picard examines aligned pairs.
# Usage: test/optical/run_optical_tests.sh <falcoAmadeus binary>
set -u
FA=${1:-./falcoAmadeus}
HERE=$(cd "$(dirname "$0")" && pwd)
W=$(mktemp -d)
trap 'rm -rf "$W"' EXIT
fail=0
python3 -I "$HERE/make_optical_data.py" 7 "$W" 1200 > /dev/null
val() { awk -F'\t' -v k="$1" '$1 == k {print $2}' "$2"; }
while IFS=$'\t' read -r D dups opt; do
  case $D in \#*) continue ;; esac
  if [ -n "${GATK:-}" ]; then
    [ -s "$W/aln.bam" ] || samtools sort -o "$W/aln.bam" "$W/aln.sam" 2> /dev/null
    "$GATK" MarkDuplicates -I "$W/aln.bam" -O "$W/md.bam" -M "$W/md.txt" \
      --OPTICAL_DUPLICATE_PIXEL_DISTANCE "$D" --VALIDATION_STRINGENCY SILENT \
      > /dev/null 2>&1
    run=$(awk -F'\t' 'p {print $7 "\t" $8; exit} /^LIBRARY/ {p = 1}' "$W/md.txt")
    [ "$run" = "$dups	$opt" ] || { echo "FAIL Picard run differs from the stored counts at $D: $run"; fail=1; }
  fi
  for t in 1 16; do
    rm -rf "$W/o"
    "$FA" -q -t $t -o "$W/o" "$W/r1.fq" --paired-with "$W/r2.fq" --optical \
      --optical-distance "$D" > /dev/null 2>&1
    R="$W/o/paired_duplication.txt"
    pairs=$(val "Total read pairs" "$R")
    groups=$(val "Groups of copies" "$R")
    got=$(val "Optical duplicates" "$R")
    if [ "$((pairs - groups))" = "$dups" ] && [ "$got" = "$opt" ]; then
      echo "ok   $D pixels t$t: $got optical of $dups duplicates, as Picard"
    else
      echo "FAIL $D pixels t$t: $got optical of $((pairs - groups)), Picard $opt of $dups"
      fail=1
    fi
  done
done < "$HERE/expected_picard.tsv"
# pairs that would not align, placed next to each other on one tile: 60
# adapter dimers, 60 pairs with a 20-base insert, 60 with mostly N in read 2;
# they are left out, so the optical count stays Picard's
python3 -I - "$W" <<'PY'
import random, sys
d = sys.argv[1]
rnd = random.Random(5)
A1 = 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTGCTTG'
A2 = 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATTAAAAAA'
rc = lambda s: s[::-1].translate(str.maketrans('ACGT', 'TGCA'))
ins = ''.join(rnd.choice('ACGT') for _ in range(20))
other = ''.join(rnd.choice('ACGT') for _ in range(100))
kinds = [('dimer', 'G' + A1, 'G' + A2), ('short', ins + A1, rc(ins) + A2),
         ('nrich', other, 'N' * 80 + other[:20])]
with open(d + '/x1.fq', 'w') as f1, open(d + '/x2.fq', 'w') as f2:
    for kind, s1, s2 in kinds:
        s1 = (s1 + 'A' * 100)[:100]
        s2 = (s2 + 'A' * 100)[:100]
        for k in range(60):
            name = 'A01234:77:HABCDEFXX:1:1101:%d:%d' % (30000 + 7 * k, 30000 + 5 * k)
            f1.write('@%s 1:N:0:ACGT\n%s\n+\n%s\n' % (name, s1, 'I' * 100))
            f2.write('@%s 2:N:0:ACGT\n%s\n+\n%s\n' % (name, s2, 'I' * 100))
PY
cat "$W/r1.fq" "$W/x1.fq" > "$W/y1.fq"
cat "$W/r2.fq" "$W/x2.fq" > "$W/y2.fq"
while IFS=$'\t' read -r D dups opt; do
  case $D in \#*) continue ;; esac
  rm -rf "$W/o"
  "$FA" -q -t 4 -o "$W/o" "$W/y1.fq" --paired-with "$W/y2.fq" --optical \
    --optical-distance "$D" > /dev/null 2>&1
  R="$W/o/paired_duplication.txt"
  got=$(val "Optical duplicates" "$R")
  out=$(awk -F'\t' '$1 ~ /^Pairs left out/ {print $2}' "$R")
  if [ "$got" = "$opt" ] && [ "$out" = 180 ]; then
    echo "ok   $D pixels: 180 pairs that would not align left out, $got optical as Picard"
  else echo "FAIL $D pixels with pairs that would not align: $got optical, $out left out"; fail=1; fi
done < "$HERE/expected_picard.tsv"
# names without coordinates
awk 'NR % 4 == 1 {$0 = "@r" NR} 1' "$W/r1.fq" > "$W/n1.fq"
awk 'NR % 4 == 1 {$0 = "@r" NR} 1' "$W/r2.fq" > "$W/n2.fq"
rm -rf "$W/o"
"$FA" -q -o "$W/o" "$W/n1.fq" --paired-with "$W/n2.fq" --optical > /dev/null 2>&1
if [ "$(val "Optical duplicates" "$W/o/paired_duplication.txt")" = 0 ] &&
   [ "$(val "Pairs with a location in the name" "$W/o/paired_duplication.txt")" = 0 ]; then
  echo "ok   names without coordinates: none located, no optical duplicates"
else echo "FAIL names without coordinates"; fail=1; fi
rm -rf "$W/o"
"$FA" -q -o "$W/o" "$W/r1.fq" --optical > /dev/null 2>&1
if [ -n "$(val "Optical duplicates" "$W/o/aware_duplication.txt")" ]; then
  echo "ok   single-end: reported in aware_duplication.txt"
else echo "FAIL single-end report"; fail=1; fi
exit $fail
