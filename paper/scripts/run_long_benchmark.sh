#!/usr/bin/env bash
# Long runs (minutes per tool) of falcoAmadeus and falco 2 on large inputs:
# wall, CPU and peak memory per run, and the resident memory sampled every
# 0.5 s, to see throughput and memory over a whole run rather than a few
# seconds.
#
# Usage: run_long_benchmark.sh <falcoAmadeus> <falco2> <out_dir> <cases.tsv>
#                              [repetitions]
# cases.tsv: one case per line, tab-separated:
#   name  threads  R1  [R2]
# (R2 given: paired; falcoAmadeus gets --paired-with R2, falco 2 both files).
# Writes <out_dir>/results.tsv and <out_dir>/rss/<case>_<tool>_t<T>_r<k>.tsv
# (seconds, resident MB); prints a summary at the end.
set -euo pipefail
FA=$1; F2=$2; OUT=$3; CASES=$4; REPS=${5:-2}
mkdir -p "$OUT/rss"
RES=$OUT/results.tsv
[ -s "$RES" ] || printf 'case\ttool\tthreads\trep\twall_s\tuser_s\tsys_s\tmax_rss_mb\tmedian_rss_mb\tinput_gb\tgb_per_s\n' > "$RES"
W=$(mktemp -d -p "$OUT")
trap 'rm -rf "$W"' EXIT

# runs "$@" under /usr/bin/time, sampling the child's resident memory
run_one() {
  local series=$1; shift
  /usr/bin/time -f "%e %U %S %M" -o "$W/tm" "$@" > /dev/null 2> "$W/err" &
  local tp=$! child="" t0 now
  t0=$(date +%s.%N)
  : > "$series"
  while kill -0 "$tp" 2> /dev/null; do
    [ -z "$child" ] && child=$(pgrep -P "$tp" | head -1 || true)
    if [ -n "$child" ] && [ -r "/proc/$child/status" ]; then
      now=$(date +%s.%N)
      awk -v t="$(echo "$now - $t0" | bc)" '/^VmRSS:/ {printf "%.1f\t%.1f\n", t, $2 / 1024}' \
        "/proc/$child/status" >> "$series" 2> /dev/null || true
    fi
    sleep 0.5
  done
  wait "$tp" || { echo "FAILED: $*" >&2; cat "$W/err" >&2; return 1; }
}

while IFS=$'\t' read -r name t r1 r2; do
  [ -z "$name" ] || [ "${name:0:1}" = "#" ] && continue
  gb=$(du -bc "$r1" ${r2:+"$r2"} | tail -1 | awk '{printf "%.2f", $1 / 1e9}')
  for rep in $(seq "$REPS"); do
    for tool in fa f2; do
      rm -rf "$W/o"
      if [ $tool = fa ]; then
        cmd=("$FA" -q -t "$t" -o "$W/o" "$r1")
        [ -n "$r2" ] && cmd+=(--paired-with "$r2")
      else
        cmd=("$F2" -t "$t" -o "$W/o" "$r1")
        [ -n "$r2" ] && cmd+=("$r2")
      fi
      s="$OUT/rss/${name}_${tool}_t${t}_r${rep}.tsv"
      run_one "$s" "${cmd[@]}"
      read -r wall user sys rss < "$W/tm"
      med=$(cut -f2 "$s" | sort -n | awk '{a[NR] = $1} END {print (NR ? a[int((NR + 1) / 2)] : 0)}')
      printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%.0f\t%.0f\t%s\t%.3f\n' "$name" $tool \
        "$t" "$rep" "$wall" "$user" "$sys" "$(echo "$rss / 1024" | bc -l)" \
        "$med" "$gb" "$(echo "$gb / $wall" | bc -l)" >> "$RES"
      echo "$name $tool t$t rep$rep: ${wall}s" >&2
    done
  done
done < "$CASES"

# summary: median over repetitions per case and tool
python3 -I - "$RES" <<'EOF'
import csv, statistics as st, sys
from collections import defaultdict
d = defaultdict(list)
order = []
for r in csv.DictReader(open(sys.argv[1]), delimiter="\t"):
    k = (r["case"], r["threads"])
    if k not in order:
        order.append(k)
    d[k + (r["tool"],)].append(r)
print(f"{'case':22s} {'t':>3s} | {'falcoAmadeus':>26s} | {'falco 2':>26s} | speed-up")
for c, t in order:
    row = []
    for tool in ("fa", "f2"):
        rs = d[(c, t, tool)]
        if not rs:
            row.append(None)
            continue
        w = st.median(float(x["wall_s"]) for x in rs)
        m = max(float(x["max_rss_mb"]) for x in rs)
        g = st.median(float(x["gb_per_s"]) for x in rs)
        row.append((w, m, g))
    fa, f2 = row
    def f(x):
        return f"{x[0]:7.1f}s {x[2]:5.2f}GB/s {x[1]:6.0f}MB" if x else " " * 26
    sp = f"{f2[0] / fa[0]:.2f}x" if fa and f2 else ""
    print(f"{c:22s} {t:>3s} | {f(fa)} | {f(f2)} | {sp}")
EOF
