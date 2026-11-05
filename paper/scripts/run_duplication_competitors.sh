#!/usr/bin/env bash
# Duplicate rates reported by alignment-free tools on the simulated library
# (simulate_duplicates.py), single-end (read 1) and paired-end, with the
# known PCR duplicate rate in truth.tsv.
#
# Usage: run_duplication_competitors.sh <dir with sim_R1/R2.fastq, truth.tsv>
#                                        <falcoAmadeus> <falco 1> <falco 2>
#                                        <fastqc> <fastp> [threads]
# Needs clumpify.sh (BBTools) and picard on PATH. Writes <dir>/results.tsv
# (tool, mode, setting, reported duplicate rate in %, wall seconds).
set -euo pipefail
D=$1; FA=$2; F1=$3; F2=$4; FQC=$5; FASTP=$6; T=${7:-16}
R1=$D/sim_R1.fastq; R2=$D/sim_R2.fastq
W=$D/work; rm -rf "$W"; mkdir -p "$W"
OUT=$D/results.tsv
printf 'tool\tmode\tsetting\tduplicate_rate_pct\tseconds\n' > "$OUT"
row() { printf '%s\t%s\t%s\t%s\t%s\n' "$1" "$2" "$3" "$4" "$5" >> "$OUT"; }
now() { date +%s.%N; }
el() { echo "$(now) - $1" | bc; }
# 100 - "Total Deduplicated Percentage" of a fastqc_data.txt
fqc_dup() { awk -F'\t' '/^#Total Deduplicated Percentage/ {printf "%.4f", 100 - $2}' "$1"; }

# ---- FastQC-type tools: one file at a time ----
for m in 1 2; do
  r=$R1; [ $m = 2 ] && r=$R2
  t0=$(now); "$FQC" -q -t "$T" --extract -o "$W" "$r" > /dev/null 2>&1
  row FastQC "SE read $m" default "$(fqc_dup "$W"/sim_R${m}_fastqc/fastqc_data.txt)" "$(el $t0)"
  t0=$(now); "$F1" -q -o "$W/f1_$m" "$r" > /dev/null 2>&1
  row "falco 1.3.2" "SE read $m" default "$(fqc_dup "$W/f1_$m/fastqc_data.txt")" "$(el $t0)"
  t0=$(now); "$F2" -t "$T" -o "$W/f2_$m" "$r" > /dev/null 2>&1
  f=$(find "$W/f2_$m" -name fastqc_data.txt | head -1)
  row "falco 2" "SE read $m" default "$(fqc_dup "$f")" "$(el $t0)"
done

# ---- fastp: exact, single-end and paired, with and without UMI ----
fastp_dup() { python3 -I -c 'import json,sys; print("%.4f" % (100 * json.load(open(sys.argv[1]))["duplication"]["rate"]))' "$1"; }
for acc in 3 6; do
  t0=$(now); "$FASTP" -w "$T" --dup_calc_accuracy $acc -i "$R1" -o /dev/null -j "$W/fp.json" -h /dev/null 2> /dev/null
  row fastp SE "accuracy $acc" "$(fastp_dup "$W/fp.json")" "$(el $t0)"
  t0=$(now); "$FASTP" -w "$T" --dup_calc_accuracy $acc -i "$R1" -I "$R2" -o "$W/fp1.fq" -O "$W/fp2.fq" -j "$W/fp.json" -h /dev/null 2> /dev/null
  row fastp PE "accuracy $acc" "$(fastp_dup "$W/fp.json")" "$(el $t0)"
  t0=$(now); "$FASTP" -w "$T" --dup_calc_accuracy $acc -U --umi_loc read1 --umi_len 8 -i "$R1" -o /dev/null -j "$W/fp.json" -h /dev/null 2> /dev/null
  row fastp SE "accuracy $acc, UMI read1 8" "$(fastp_dup "$W/fp.json")" "$(el $t0)"
  t0=$(now); "$FASTP" -w "$T" --dup_calc_accuracy $acc -U --umi_loc per_read --umi_len 8 -i "$R1" -I "$R2" -o "$W/fp1.fq" -O "$W/fp2.fq" -j "$W/fp.json" -h /dev/null 2> /dev/null
  row fastp PE "accuracy $acc, UMI per read 8" "$(fastp_dup "$W/fp.json")" "$(el $t0)"
done
rm -f "$W/fp1.fq" "$W/fp2.fq"

# ---- clumpify: dedupe with substitutions, removed reads / input reads ----
clump_dup() { awk '/^Reads In:/ {i=$3} /^Duplicates Found:/ {d=$3} END {printf "%.4f", 100*d/i}' "$1"; }
for s in 0 2; do
  t0=$(now); clumpify.sh -Xmx16g threads="$T" in="$R1" out="$W/c1.fq" dedupe subs=$s overwrite=t > "$W/c.log" 2>&1
  row clumpify SE "subs=$s" "$(clump_dup "$W/c.log")" "$(el $t0)"
  t0=$(now); clumpify.sh -Xmx16g threads="$T" in="$R1" in2="$R2" out="$W/c1.fq" out2="$W/c2.fq" dedupe subs=$s overwrite=t > "$W/c.log" 2>&1
  row clumpify PE "subs=$s" "$(clump_dup "$W/c.log")" "$(el $t0)"
done
rm -f "$W/c1.fq" "$W/c2.fq"

# ---- Picard EstimateLibraryComplexity: unaligned pairs, mismatch rate,
# optionally the UMI as barcode tag (UMI bases removed from read 1) ----
python3 -I - "$R1" "$R2" "$W" <<'EOF'
import sys
r1, r2, w = sys.argv[1:4]
hdr = "@HD\tVN:1.6\tSO:queryname\n@RG\tID:A\tSM:sim\tLB:lib\tPL:ILLUMINA\n"
with open(r1) as a, open(r2) as b, open(w + "/pe.sam", "w") as o, \
        open(w + "/pe_umi.sam", "w") as u, open(w + "/se.sam", "w") as s:
    for f in (o, u, s):
        f.write(hdr)
    while True:
        h = a.readline()
        if not h:
            break
        s1 = a.readline().strip(); a.readline(); q1 = a.readline().strip()
        b.readline(); s2 = b.readline().strip(); b.readline(); q2 = b.readline().strip()
        n = h[1:].strip()
        o.write(f"{n}\t77\t*\t0\t0\t*\t*\t0\t0\t{s1}\t{q1}\tRG:Z:A\n")
        o.write(f"{n}\t141\t*\t0\t0\t*\t*\t0\t0\t{s2}\t{q2}\tRG:Z:A\n")
        u.write(f"{n}\t77\t*\t0\t0\t*\t*\t0\t0\t{s1[8:]}\t{q1[8:]}\tRG:Z:A\tBX:Z:{s1[:8]}\n")
        u.write(f"{n}\t141\t*\t0\t0\t*\t*\t0\t0\t{s2[8:]}\t{q2[8:]}\tRG:Z:A\tBX:Z:{s2[:8]}\n")
        s.write(f"{n}\t4\t*\t0\t0\t*\t*\t0\t0\t{s1}\t{q1}\tRG:Z:A\n")
EOF
elc_dup() { awk -F'\t' 'f && NF > 5 {printf "%.4f", 100 * $(c); exit} /^LIBRARY/ {for (i = 1; i <= NF; i++) if ($i == "PERCENT_DUPLICATION") c = i; f = 1}' "$1"; }
t0=$(now); picard -Xmx16g EstimateLibraryComplexity I="$W/pe.sam" O="$W/elc.txt" > "$W/elc.log" 2>&1
row "Picard ELC" PE "MAX_DIFF_RATE 0.03" "$(elc_dup "$W/elc.txt")" "$(el $t0)"
t0=$(now); picard -Xmx16g EstimateLibraryComplexity I="$W/pe_umi.sam" O="$W/elc_umi.txt" READ_ONE_BARCODE_TAG=BX READ_TWO_BARCODE_TAG=BX > "$W/elc_umi.log" 2>&1
row "Picard ELC" PE "MAX_DIFF_RATE 0.03, UMIs as read barcode tags" "$(elc_dup "$W/elc_umi.txt")" "$(el $t0)"
t0=$(now); if picard -Xmx16g EstimateLibraryComplexity I="$W/se.sam" O="$W/elc_se.txt" > "$W/elc_se.log" 2>&1; then
  row "Picard ELC" SE "MAX_DIFF_RATE 0.03" "$(elc_dup "$W/elc_se.txt")" "$(el $t0)"
else
  row "Picard ELC" SE "MAX_DIFF_RATE 0.03" "failed" "$(el $t0)"
fi

# ---- falcoAmadeus: exact, two mismatches, and with the UMI ----
fa_val() { awk -F'\t' -v k="$2" '$1 == k {printf "%.4f", $2}' "$1"; }
for set in "default (FastQC estimate; PE: exact pairs)||" "2 mismatches|--dup-mismatches 2|" "2 mismatches, UMIs 1:8 (1 mismatch)|--dup-mismatches 2 --umi-r1 1:8 --umi-mismatches 1|--umi-r2 1:8"; do
  name=${set%%|*}; rest=${set#*|}; args=${rest%%|*}; pe_extra=${rest#*|}
  rm -rf "$W/fa"
  t0=$(now); "$FA" -q -t "$T" -o "$W/fa" "$R1" $args > /dev/null 2>&1
  if [ -z "$args" ]; then v=$(fqc_dup "$W/fa/fastqc_data.txt"); else v=$(fa_val "$W/fa/aware_duplication.txt" Aware); fi
  row falcoAmadeus SE "$name" "$v" "$(el $t0)"
  rm -rf "$W/fa"
  t0=$(now); "$FA" -q -t "$T" -o "$W/fa" "$R1" --paired-with "$R2" $args $pe_extra > /dev/null 2>&1
  if [ -z "$args" ]; then v=$(fa_val "$W/fa/paired_duplication.txt" "Paired-end exact"); else v=$(fa_val "$W/fa/paired_duplication.txt" "Paired-end aware"); fi
  row falcoAmadeus PE "$name" "$v" "$(el $t0)"
done
column -t -s$'\t' "$OUT"
cat "$D/truth.tsv"
