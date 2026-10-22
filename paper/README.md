# Manuscript figures and measurements

Scripts and measured data behind the figures and supporting tables of the
falcoAmadeus manuscript.

- `scripts/make_figures.py <data> <out>`: one PDF per figure panel and the
  table values (`table_*.tsv`).
- `scripts/make_supporting_tables.py <data> <panels> <out>`: S1 to S3 Tables.
- `scripts/run_duplication_simulation.sh <falcoAmadeus> <out.tsv>`: paired
  duplicate rates on simulated libraries (`data/duplication_simulation.tsv`).

Data (`data/`), measured on an AMD Ryzen Threadripper 7970X:

- `se_falcoamadeus.tsv`, `se_falco2.tsv`, `se_parent_tools.tsv`: single-end
  benchmark (wall time, CPU time, peak memory, exit code, reads counted).
- `pe_competitors_and_twofile.tsv`, `pe_falcoamadeus_paired.tsv`: paired-end
  benchmark on ten million HG002 pairs (SRR12898355).
- `assay_validation_runs.tsv`, `assay_validation_grades.tsv`: the 120 public
  runs and their module grades (validation/assay_profiles).
