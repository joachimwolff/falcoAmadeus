#!/usr/bin/env bash
# Machine-readable summary (falcoAmadeus_summary.json, _mqc.json): valid
# JSON for a single-end run with --complexity and a paired run; the Q20/Q30
# base shares equal a direct count over the quality strings; the paired
# summary holds both mates, the paired duplication, complexity and insert
# size; the MultiQC file is general-statistics custom content.
# Usage: test/json/run_json_tests.sh <falcoAmadeus binary>
set -u
FA=${1:-./falcoAmadeus}
W=$(mktemp -d)
trap 'rm -rf "$W"' EXIT
fail=0
# 5000 pairs, 100 bases, qualities cycling through Phred 2..41
awk -v o1="$W/r1.fq" -v o2="$W/r2.fq" '
function rnd() { x = (x * 16807) % 2147483647; return x / 2147483647 }
BEGIN {
  x = 5
  for (i = 0; i < 5000; i++) {
    s1 = ""; s2 = ""; q = ""
    for (k = 0; k < 100; k++) {
      s1 = s1 substr("ACGT", int(rnd() * 4) + 1, 1)
      s2 = s2 substr("ACGT", int(rnd() * 4) + 1, 1)
      q = q sprintf("%c", 35 + (i + k) % 40)
    }
    if (i > 0 && rnd() < 0.4) { s1 = p1; s2 = p2 }  # copies, geometric
    p1 = s1; p2 = s2
    printf "@p%d/1\n%s\n+\n%s\n", i, s1, q > o1
    printf "@p%d/2\n%s\n+\n%s\n", i, s2, q > o2
  }
}'
"$FA" -q --complexity -o "$W/se" "$W/r1.fq" > /dev/null 2>&1
"$FA" -q -o "$W/pe" "$W/r1.fq" --paired-with "$W/r2.fq" > /dev/null 2>&1
python3 -I - "$W" <<'PY' || fail=1
import json, sys
w = sys.argv[1]
bad = 0
def ok(c, msg):
    global bad
    print(("ok   " if c else "FAIL ") + msg)
    bad |= not c
# direct count over the qualities of read 1
q20 = q30 = n = 0
for i, line in enumerate(open(w + "/r1.fq")):
    if i % 4 == 3:
        for c in line.strip():
            p = ord(c) - 33
            n += 1; q20 += p >= 20; q30 += p >= 30
se = json.load(open(w + "/se/falcoAmadeus_summary.json"))
v = se["views"]["reads"]
ok(abs(v["q30_bases_pct"] - 100 * q30 / n) < 1e-4 and
   abs(v["q20_bases_pct"] - 100 * q20 / n) < 1e-4, "single-end Q20/Q30")
ok(v["basic_statistics"]["total_sequences"] == 5000, "single-end total sequences")
c = se["complexity"]  # the curve, or why it cannot be extrapolated
ok("10x" in c.get("complexity", {}) or "extrapolation" in c,
   "single-end complexity")
pe = json.load(open(w + "/pe/falcoAmadeus_summary.json"))
ok(set(pe["views"]) >= {"forward", "reverse"}, "paired views")
ok(pe["paired_duplication"]["total_read_pairs"] == 5000, "paired duplication")
ok("insert_size_peak" in pe["insert_size"], "insert size")
m = json.load(open(w + "/pe/falcoAmadeus_mqc.json"))
ok(m["plot_type"] == "generalstats" and "q30_bases_pct" in
   list(m["data"].values())[0], "MultiQC general statistics")
sys.exit(bad)
PY
exit $fail
