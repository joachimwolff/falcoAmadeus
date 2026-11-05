#!/usr/bin/env bash
# Assay-profile validation: for every public run in runs.tsv, grade the clean
# reads and each injected defect under the generic profile and under the
# run's assay profile, and record every module grade.
#
# Usage: run_validation.sh <falcoAmadeus> <runs.tsv> <fastq_dir> <out.tsv>
#   fastq_dir holds <run>.fastq (first 1M reads of R1 for paired runs)
set -u
FA=$1; RUNS=$2; FQ=$3; OUT=$4
HERE=$(cd "$(dirname "$0")" && pwd)
W=$(mktemp -d)
trap 'rm -rf "$W"' EXIT
DEFECTS="adapter qualtail midbias gccontam dupspike ncontent"
FRACTION=0.5

profile_of() {  # runs.tsv assay label -> --assay value
  case $1 in
    amplicon16s) echo 16s ;;
    *) echo "$1" ;;
  esac
}

echo -e "assay\trun\tcondition\tprofile\tmodule\tgrade" > "$OUT"
grade() {  # assay run condition profile fastq
  local o="$W/o"
  rm -rf "$o"; mkdir -p "$o"
  "$FA" -t 4 -q --assay "$4" -o "$o" "$5" > /dev/null 2>&1 || {
    echo -e "$1\t$2\t$3\t$4\tRUN_FAILED\tNA" >> "$OUT"; return; }
  awk -F'\t' -v a="$1" -v r="$2" -v c="$3" -v p="$4" \
    '{print a"\t"r"\t"c"\t"p"\t"$2"\t"$1}' "$o/summary.txt" >> "$OUT"
}

tail -n +2 "$RUNS" | while IFS=$'\t' read -r assay run _rest; do
  f="$FQ/$run.fastq"
  [ -s "$f" ] || { echo "missing $f" >&2; continue; }
  prof=$(profile_of "$assay")
  for p in generic "$prof"; do
    [ "$p" = generic ] && [ "$prof" = generic ] && continue
    grade "$assay" "$run" clean "$p" "$f"
  done
  for d in $DEFECTS; do
    python3 -I "$HERE/inject_defect.py" "$f" "$W/d.fastq" "$d" "$FRACTION" --seed 7
    for p in generic "$prof"; do
      grade "$assay" "$run" "$d" "$p" "$W/d.fastq"
    done
  done
  echo "done $assay $run" >&2
done
