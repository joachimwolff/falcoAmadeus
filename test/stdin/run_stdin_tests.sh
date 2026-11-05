#!/usr/bin/env bash
# FASTQ from a pipe: the reports equal those of the same file read directly
# (apart from the file name) at one and several threads, also for a file
# whose last line has no newline.
# Usage: test/stdin/run_stdin_tests.sh <falcoAmadeus binary>
set -u
FA=${1:-./falcoAmadeus}
W=$(mktemp -d)
trap 'rm -rf "$W"' EXIT
fail=0

# 30000 reads of 76 bases, a few of them repeated
awk 'BEGIN {
  q = "#"; for (p = 2; p <= 76; p++) q = q "I"
  for (i = 0; i < 30000; i++) {
    s = ""
    for (p = 1; p <= 76; p++)
      s = s substr("ACGT", ((i % 997) * 7 + p * 13 + int(p / 5)) % 4 + 1, 1)
    printf "@r%d\n%s\n+\n%s\n", i, s, q
  }
}' > "$W/reads.fastq"
head -c -1 "$W/reads.fastq" > "$W/nonl.fastq"  # no newline at the end

data() { grep -v '^Filename' "$1/fastqc_data.txt" | tail -n +2; }

for f in reads nonl; do
  rm -rf "$W/file"; mkdir "$W/file"
  "$FA" -q -o "$W/file" "$W/$f.fastq" > /dev/null 2>&1
  for t in 1 8; do
    rm -rf "$W/pipe"; mkdir "$W/pipe"
    cat "$W/$f.fastq" | "$FA" -q -t "$t" -f fastq -o "$W/pipe" - \
      > /dev/null 2>&1
    if diff <(data "$W/file") <(data "$W/pipe") > /dev/null &&
       diff <(cut -f1,2 "$W/file/summary.txt") \
            <(cut -f1,2 "$W/pipe/summary.txt") > /dev/null; then
      echo "ok   $f -t $t: pipe = file"
    else
      echo "FAIL $f -t $t: pipe differs from file"; fail=1
    fi
  done
done

exit $fail
