# Manuscript figures and measurements

Scripts and measured data behind the figures and supporting tables of the
falcoAmadeus manuscript.

- `scripts/make_figures.py <data> <out>`: one PDF per figure panel and the
  table values (`table_*.tsv`).
- `scripts/make_supporting_tables.py <data> <panels> <out>`: S1 to S5 Tables.
- `scripts/run_duplication_simulation.sh <falcoAmadeus> <out.tsv>`: paired
  duplicate rates on simulated libraries (`data/duplication_simulation.tsv`).

Data (`data/`), measured on an AMD Ryzen Threadripper 7970X:

- `se_falcoamadeus.tsv`, `se_falco2.tsv`, `se_parent_tools.tsv`: single-end
  benchmark (wall time, CPU time, peak memory, exit code, reads counted).
- `pe_competitors_and_twofile.tsv`, `pe_falcoamadeus_paired.tsv`: paired-end
  benchmark on ten million HG002 pairs (SRR12898355).
- `assay_validation_runs.tsv`, `assay_validation_grades.tsv`: the 120 public
  runs and their module grades (validation/assay_profiles).
- `duplication_competitors.tsv`, `duplication_competitors_truth.tsv`:
  duplicate rates of FastQC 0.12.1, falco 1.3.2, falco 2.0.3, fastp 0.24.3,
  clumpify (BBTools 39.01), Picard EstimateLibraryComplexity 2.20.4 and
  falcoAmadeus on a simulated library with known duplication
  (`scripts/simulate_duplicates.py`, `scripts/run_duplication_competitors.sh`).
- `picard_comparison_<set>.tsv`, `picard_groups_<set>.tsv`,
  `picard_set_mismatches_<set>.tsv`: falcoAmadeus against Picard
  MarkDuplicates on the same pairs, exact (`scripts/compare_picard.py`) and
  with 1 and 2 mismatches (`scripts/compare_picard_groups.py` on the groups
  written with `FA_AWARE_CLUSTERS`).
- `optical_validation/`, `trim_preview_validation/`,
  `complexity_insert_validation/`, `screen_validation/`,
  `umi_qc_validation/`: comparisons with Picard, cutadapt, preseq, fastp and
  an independent UMI computation, each with its README.
