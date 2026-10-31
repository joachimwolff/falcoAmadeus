#!/usr/bin/env bash
# Contamination screen: builds the 15-reference panel from the downloaded
# references (references_manifest.tsv: name, file, ... in <refs_dir>),
# simulates 1,000,000 reads in known proportions (simulate_mixture.py; human
# and mouse reads from chromosome 1, the database covers whole genomes) and
# screens them, then screens public runs given as name=fastq.
#
# Usage: run_screen_validation.sh <falcoAmadeus> <refs_dir> <out_dir>
#                                 [name=fastq ...]
set -euo pipefail
FA=$1; REFS=$2; OUT=$3; shift 3
HERE=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$OUT"
R=""
while IFS=$'\t' read -r name file rest; do
  [ "$name" = name ] && continue
  R="$R,$name=$REFS/$file"
done < "$REFS/MANIFEST.tsv"
"$FA" -t 16 --screen-build "$OUT/panel.fascreen" --screen-refs "${R#,}"
zcat "$REFS"/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz |
  awk '/^>/ {p = ($1 == ">1")} p' > "$OUT/hs_chr1.fa"
zcat "$REFS"/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz |
  awk '/^>/ {p = ($1 == ">1")} p' > "$OUT/mm_chr1.fa"
python3 -I "$HERE/simulate_mixture.py" 1000000 3 "$OUT/mix.fastq" \
  human=0.90="$OUT/hs_chr1.fa" mouse=0.04="$OUT/mm_chr1.fa" \
  ecoli=0.02="$REFS/ecoli.fa.gz" phix=0.01="$REFS/phix.fa.gz" \
  mycoplasma=0.01="$REFS/mycoplasma.fa.gz" \
  human_rrna=0.01="$REFS/rrna_human.fa.gz" \
  vectors=0.005="$REFS/univec_core.fa.gz" lambda=0.005="$REFS/lambda.fa.gz"
"$FA" -q -t 16 --screen "$OUT/panel.fascreen" -o "$OUT/mix" "$OUT/mix.fastq"
for nf in "$@"; do
  "$FA" -q -t 16 --screen "$OUT/panel.fascreen" -o "$OUT/run_${nf%%=*}" \
    "${nf#*=}"
done
rm -f "$OUT/hs_chr1.fa" "$OUT/mm_chr1.fa" "$OUT/mix.fastq"
