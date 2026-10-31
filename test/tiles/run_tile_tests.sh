#!/usr/bin/env bash
# Per tile sequence quality when read names change format after the first
# read (make_mixed_names.py): the tile field then holds arbitrary numbers,
# some above 32 bits. The section equals stock falco 1.3.1's
# (expected_per_tile.tsv) for plain and gzip input, 1 and 4 threads, and the
# forward mate of a pair.
# Usage: test/tiles/run_tile_tests.sh <falcoAmadeus binary>
set -u
FA=${1:-./falcoAmadeus}
HERE=$(cd "$(dirname "$0")" && pwd)
W=$(mktemp -d)
trap 'rm -rf "$W"' EXIT
fail=0
python3 -I "$HERE/make_mixed_names.py" "$W/m.fq"
gzip -1 -c "$W/m.fq" > "$W/m.fq.gz"
cp "$W/m.fq" "$W/m2.fq"
section() { awk '/>>Per tile/,/>>END/' "$1" 2> /dev/null; }
for input in m.fq m.fq.gz; do
  for t in 1 4; do
    rm -rf "$W/o"
    "$FA" -q -t $t -o "$W/o" "$W/$input" > /dev/null 2>&1
    rc=$?
    if [ $rc = 0 ] && cmp -s <(section "$W/o/fastqc_data.txt") "$HERE/expected_per_tile.tsv"; then
      echo "ok   $input t$t: per tile quality equals falco"
    else echo "FAIL $input t$t (exit $rc)"; fail=1; fi
  done
done
rm -rf "$W/p"
"$FA" -q -t 4 -o "$W/p" "$W/m.fq" --paired-with "$W/m2.fq" > /dev/null 2>&1
if cmp -s <(section "$W/p/m_forward/fastqc_data.txt") "$HERE/expected_per_tile.tsv"; then
  echo "ok   paired forward mate: per tile quality equals falco"
else echo "FAIL paired forward mate"; fail=1; fi
exit $fail
