#!/usr/bin/env bash
# Optical duplicates (--optical) on read pairs from make_optical_data.py:
# copies are sequence-identical and at the same coordinates, so the groups
# of copies equal Picard's duplicate sets. Duplicates and optical duplicates
# equal those of Picard MarkDuplicates (GATK 4.6.2.0) on the matching
# alignments at 50, 100 and 2500 pixels (expected_picard.tsv; with GATK=gatk
# and samtools on PATH also run directly), for 1 and 16 threads. Names
# without coordinates give no optical duplicates, and single-end input is
# reported per read.
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
