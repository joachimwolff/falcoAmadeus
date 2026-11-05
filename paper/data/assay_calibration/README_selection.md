# Run selection for QC-tool calibration (assay_calibration)

Date of search: 2026-10-30. Calibration set, disjoint from the validation set in `../assay_validation/runs_final.tsv`.

Output: `runs.tsv` (70 rows: 8 runs for each of rnaseq, smallrna, atac, chipseq, wgbs, rrbs, amplicon16s, exome; 6 runs for emseq). Each FASTQ is the first 1,000,000 spots of read 1 (`fastq/<run>.fastq`, 4,000,000 lines each). `download_report.tsv` lists the line count of every FASTQ and the one replacement.

## Sources and queries

NCBI E-utilities, db=sra, followed by a runinfo check of every accepted run (`https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/runinfo?acc=<RUN>`). Candidate records were parsed with esummary (experiment title, study name, organism, platform, LibraryStrategy, LibraryLayout, BioProject, BioSample, instrument, per-run total_spots and is_public). Searches used the same LibraryStrategy and platform terms as the validation README, with one esearch per release year (2012 to 2026, up to 150 IDs per year and query), or without a date window where the pool is small.

- rnaseq: `"RNA-Seq"[Strategy] AND "TRANSCRIPTOMIC"[Source] AND "Homo sapiens"[Organism] AND "ILLUMINA"[Platform]`
- smallrna: `"miRNA-Seq"[Strategy] AND "Homo sapiens"[Organism] AND "ILLUMINA"[Platform]`
- atac: `"ATAC-seq"[Strategy] AND "Homo sapiens"[Organism] AND "ILLUMINA"[Platform]`
- chipseq: `"ChIP-Seq"[Strategy] AND "Homo sapiens"[Organism] AND "ILLUMINA"[Platform]`
- wgbs: `"Bisulfite-Seq"[Strategy] AND "Homo sapiens"[Organism] AND ("WGBS" OR "whole genome bisulfite")[All Fields] AND "ILLUMINA"[Platform]`
- emseq: `"Bisulfite-Seq"[Strategy] AND ("EM-seq" OR "EMseq" OR "enzymatic methyl-seq" OR "enzymatic methylation" OR "enzymatic conversion")[All Fields] AND ("Homo sapiens" OR "Mus musculus")[Organism] AND "ILLUMINA"[Platform]`, plus the broader per-year queries without organism restriction (results were then restricted to human and mouse, see below)
- rrbs: `"RRBS"[All Fields] AND "ILLUMINA"[Platform]`; `"reduced representation bisulfite"[All Fields] AND "ILLUMINA"[Platform]`
- amplicon16s: `"AMPLICON"[Strategy] AND ("16S" OR "16S rRNA" OR "rDNA" OR "microbiome" OR "bacteria" OR "Targeted Locus")[All Fields] AND "ILLUMINA"[Platform]`
- exome: `"WXS"[Strategy] AND "Homo sapiens"[Organism] AND "ILLUMINA"[Platform]`

## Filters

Applied to the esummary records:

1. Platform ILLUMINA; LibraryStrategy equal to the assay target (RNA-Seq with LibrarySource TRANSCRIPTOMIC for rnaseq; miRNA-Seq; ATAC-seq; ChIP-Seq; Bisulfite-Seq for wgbs, emseq and rrbs; AMPLICON; WXS).
2. Organism: Homo sapiens for rnaseq, smallrna, atac, chipseq, wgbs and exome. For emseq and rrbs, Homo sapiens or Mus musculus only. amplicon16s is not restricted by organism (microbiome).
3. Title or study text must match the assay: EM-seq / enzymatic methyl for emseq; RRBS / reduced representation bisulfite for rrbs; WGBS / whole genome bisulfite (and no EM-seq or RRBS terms) for wgbs; for amplicon16s, a bacterial 16S marker or microbiome term (16S, rRNA, rDNA, Bv6, microbio, bacteri, gut, soil, marine, fecal, stool, metagenom) and no eukaryotic or editing term (18S, ITS, editing, CRISPR, saturation, sgRNA, guide).
4. At least one run with total_spots >= 1,000,000 and is_public = true.
5. Studies whose title is a GEO placeholder ("currently private") were ranked last and not used.

Each accepted run was then checked with runinfo: LibraryStrategy equal to the assay target, Platform ILLUMINA, spots >= 1,000,000, avgLength >= 36, empty dbgap_study_accession and Consent public (controlled-access runs were rejected before download), BioProject non-empty, ScientificName Homo sapiens where required, and LibrarySource TRANSCRIPTOMIC for rnaseq. Runs that failed a check were rejected and the next candidate was taken.

## Selection rules

- Candidates were ordered round-robin over BioProjects, most recent first, so that each BioProject contributes at most 2 runs per assay. Across all 70 rows, no BioProject appears more than twice.
- The first 8 verified runs (6 for emseq, see below) were taken in that order.

## Exclusion of the validation set

Excluded from every search: all runs and BioProjects in `../assay_validation/runs_final.tsv` (120 runs, 81 BioProjects). The earlier validation pool `../assay_validation/runs.tsv` and its replacement list `runs_replacements.tsv` were excluded in the same way, so that no run from the removed or replaced validation rows is used either. The overlap was checked after writing `runs.tsv`: no run and no BioProject is shared with `runs_final.tsv` (comparison of the run and bioproject columns).

## Shortfall: emseq has 6 runs

The EM-seq pool that meets the organism rule (Homo sapiens or Mus musculus), the title rule, the 1,000,000-spot rule and the exclusion of validation BioProjects contains only three BioProjects: PRJNA591788 (human, 2 runs), PRJNA884710 (mouse, 2 runs) and PRJNA872013 (mouse, 2 runs). With at most 2 runs per BioProject this gives 6 runs. Non-human EM-seq runs (plant, insect, fish) were not used because the organism rule allows only human and mouse. Two human or mouse EM-seq runs remain missing; the other human and mouse EM-seq BioProjects (PRJNA200694, PRJNA600820, PRJNA834955, PRJNA1446745, PRJNA1311825) are in the validation set and are excluded. The emseq assay therefore has 6 runs, not 8.

## Notes on amplicon16s

Most AMPLICON runs in SRA with at least 1,000,000 spots are labelled "Targeted Locus (Loci)" with an environmental organism. Bacterial 16S status for the PRJNA262654, PRJNA207080 and PRJNA262656 runs rests on the experiment title tags (`Bv6`, the bacterial V6 16S marker label), not on a primer field in the record. PRJNA264350 ("16S rDNA V6") and PRJNA271708 ("16S rRNA deep sequencing") carry explicit 16S text. PRJEB36282 (18S rDNA) was excluded. SRR981224 (PRJNA219042) is a human gut microbiota run (single-end, 101 bp). The microbiome amplicon pool is small, so the set is dominated by environmental samples.

## Replacement and failures

- SRR15114598 (rrbs, PRJNA556128): `fastq-dump -X 1000000` returned 989,934 spots (3,959,736 lines), below the 1,000,000 required. It was removed and replaced with the first verified reserve SRR10336785 (rrbs, PRJNA578952, Homo sapiens, 28,565,472 spots, avgLength 300). The replacement downloaded with exactly 4,000,000 lines. The failure is recorded in `download_report.tsv`.
- A first download attempt failed for all 70 runs because of a line-ending error in the job file (the run names carried a trailing newline). It was fixed and the batch was rerun from the same job list; those failures are not counted in `download_report.tsv`.

## Download

`fastq-dump -X 1000000 --split-files` from `~/miniconda3/envs/sratoolkit/bin`, up to 8 runs in parallel, each run in its own subdirectory of a temporary directory under this folder that was removed after the run. Read 1 (`<run>_1.fastq`) for paired runs; `<run>.fastq` for single-end runs, renamed to `fastq/<run>.fastq`. Nothing was written under `/mnt/processed`.
