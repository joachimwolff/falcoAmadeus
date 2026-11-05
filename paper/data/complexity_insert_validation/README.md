# Library complexity and insert size against preseq and fastp

`paper/scripts/complexity_insert_validation/run_complexity_insert.sh`
(falcoAmadeus 4, preseq 3.2.0, fastp 0.24.3) on all read pairs of two
public runs: DRR041223 (mouse RNA-seq, 3,183,182 pairs, 2x100) and the
first 10,000,000 pairs of SRR12898355 (GIAB HG002 WGS, 2x126).

- `summary.tsv`: the yield curve of falcoAmadeus (`FA_COMPLEXITY_CURVE`)
  against `preseq lc_extrap -Q -H` on the copy-number histogram of the same
  run (`FA_COPY_HIST`), all 10,000 points; the insert-size histogram of
  `paired_insert_size.txt` (`--insert-size-pairs 0`, all pairs) against
  fastp `-w 1 -A -G` (one worker, so fastp counts every pair), all 513
  bins. Both equal on both runs.
- `<run>.curve.tsv`: the yield curve (equal to preseq's output).
- `<run>.yield.tsv`: the complexity rows of `paired_duplication.txt`.
- `<run>.insert_size.tsv`: `paired_insert_size.txt` (file names shortened).
