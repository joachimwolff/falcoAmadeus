#!/usr/bin/env python3
"""Supporting tables S1 to S5 of the falcoAmadeus manuscript as Excel files.

Usage: make_supporting_tables.py <data_dir> <panels_dir> <out_dir>
  panels_dir holds table_se.tsv, table_pe.tsv, table_assay_modules.tsv and
  table_assay_detection.tsv written by make_figures.py; data_dir holds the
  Picard comparisons (picard_comparison_<set>.tsv, picard_groups_<set>.tsv,
  picard_set_mismatches_<set>.tsv) and the comparisons with other programs
  (optical_validation, trim_preview_validation, complexity_insert_validation,
  screen_validation).
"""
import csv
import os
import sys

from openpyxl import Workbook
from openpyxl.styles import Font

DATA, PANELS, OUT = sys.argv[1:4]
os.makedirs(OUT, exist_ok=True)

TOOL = {"fa": "falcoAmadeus", "fa_pe": "falcoAmadeus paired mode",
        "f2": "falco 2.0.3", "f2o": "falco 2.0.3 --orig-dups",
        "f1": "falco 1.3.2", "fastqc": "FastQC 0.12.1"}
SHAPE = {"wide": "32 cores", "narrow": "4 cores (8 threads)"}
CACHE = {"warm": "cached", "cold": "not cached"}


def tsv(path):
    with open(path) as f:
        return list(csv.DictReader(f, delimiter="\t"))


def sheet(wb, title, header, rows, first=False):
    ws = wb.active if first else wb.create_sheet()
    ws.title = title
    ws.append(header)
    for c in ws[1]:
        c.font = Font(bold=True)
    for r in rows:
        ws.append(r)
    for col in ws.columns:
        width = max(len(str(c.value)) for c in col if c.value is not None)
        ws.column_dimensions[col[0].column_letter].width = min(width + 2, 60)


def num(x):
    try:
        return float(x)
    except ValueError:
        return x


# S1 Table: benchmarks
wb = Workbook()
se = tsv(os.path.join(PANELS, "table_se.tsv"))
sheet(wb, "Single-end",
      ["Input file", "Machine", "Page cache", "Threads", "Program",
       "Wall time (s)", "Peak memory (MB)", "Completed"],
      [[r["file"], SHAPE[r["shape"]], CACHE[r["cache"]], int(r["threads"]),
        TOOL[r["tool"]], num(r["wall_s"]), num(r["rss_mb"]),
        "yes" if r["ok"] == "1" else "no (error)"]
       for r in se], first=True)
pe = tsv(os.path.join(PANELS, "table_pe.tsv"))
sheet(wb, "Paired-end",
      ["Format", "Machine", "Threads", "Program", "Wall time (s)",
       "Peak memory (MB)", "Completed"],
      [[r["format"], SHAPE[r["shape"]], int(r["threads"]), TOOL[r["tool"]],
        num(r["wall_s"]), num(r["rss_mb"]),
        "yes" if r["ok"] == "1" else "no (error)"] for r in pe])
wb.save(os.path.join(OUT, "S1_Table.xlsx"))

# S2 Table: runs of the assay validation
runs = tsv(os.path.join(DATA, "assay_validation_runs.tsv"))
wb = Workbook()
sheet(wb, "Runs",
      ["Assay", "Run", "BioProject", "BioSample", "Library strategy",
       "Library selection", "Layout", "Reads", "Average length",
       "Instrument", "Study title"],
      [[r["assay"], r["run"], r["bioproject"], r["biosample"],
        r["library_strategy"], r["library_selection"], r["layout"],
        num(r["spots"]), num(r["avg_length"]), r["instrument"],
        r["study_title"]] for r in runs], first=True)
wb.save(os.path.join(OUT, "S2_Table.xlsx"))

# S3 Table: grades of the assay validation
wb = Workbook()
mods = tsv(os.path.join(PANELS, "table_assay_modules.tsv"))
sheet(wb, "Unmodified libraries",
      ["Assay", "Module", "Runs flagged, generic", "Runs flagged, profile",
       "Runs"],
      [[r["assay"], r["module"], int(r["generic_flagged"]),
        int(r["profile_flagged"]), int(r["runs"])] for r in mods],
      first=True)
det = tsv(os.path.join(PANELS, "table_assay_detection.tsv"))
LABEL = {"adapter": "Adapter (Adapter Content)",
         "qualtail": "Quality tail (Per base sequence quality)",
         "gccontam": "GC contamination (Per sequence GC content)",
         "dupspike": "Duplicate spike (Sequence Duplication Levels)",
         "ncontent": "N block (Per base N content)",
         "midbias": "Composition bias (Per base sequence content)"}
sheet(wb, "Injected defects",
      ["Defect (targeted module)", "Runs detected, generic",
       "Runs detected, profile", "Runs"],
      [[LABEL[r["defect"]], int(r["generic_detected"]),
        int(r["profile_detected"]), int(r["runs"])]
       for r in det if r["defect"] in LABEL])
wb.save(os.path.join(OUT, "S3_Table.xlsx"))

# S4 Table: paired duplicate rate against Picard MarkDuplicates
SETS = [("hg002", "HG002 WGS (SRR12898355)"),
        ("drr041223", "Mouse RNA-seq (DRR041223)")]
CMP = [("pairs_total", "Read pairs"),
       ("pairs_compared", "Pairs compared (both mates primary, properly "
                          "paired, mapping quality >= 30)"),
       ("picard_duplicate_pairs", "Picard duplicate pairs"),
       ("picard_duplicate_rate_pct", "Picard duplicate rate (%)"),
       ("sequence_duplicate_pairs", "falcoAmadeus duplicate pairs"),
       ("sequence_duplicate_rate_pct", "falcoAmadeus duplicate rate (%)"),
       ("sequence_extra_copies_flagged_by_picard_pct",
        "falcoAmadeus duplicates that are Picard duplicates (%)"),
       ("picard_flagged_with_sequence_partner_pct",
        "Picard duplicates with a sequence-identical pair (%)")]
cmp_rows = {k: {r["metric"]: r["value"] for r in
                tsv(os.path.join(DATA, f"picard_comparison_{k}.tsv"))}
            for k, _ in SETS}
wb = Workbook()
sheet(wb, "Comparison", ["Measure"] + [label for _, label in SETS],
      [[label] + [num(cmp_rows[k][m]) for k, _ in SETS] for m, label in CMP],
      first=True)
GROUPS = [("pairs_compared", "Pairs compared"),
          ("picard_duplicate_rate_pct", "Picard duplicate rate (%)"),
          ("falcoamadeus_duplicate_rate_pct",
           "falcoAmadeus duplicate rate (%)"),
          ("falcoamadeus_extra_copies", "falcoAmadeus additional copies"),
          ("extra_copies_flagged_by_picard_pct",
           "falcoAmadeus copies that are Picard duplicates (%)"),
          ("picard_duplicates_found_pct",
           "Picard duplicates in a falcoAmadeus group (%)")]
rows = []
for k, label in SETS:
    ex = cmp_rows[k]
    found = 100.0 * float(ex["sequence_extra_copies_flagged_by_picard"]) / \
        float(ex["picard_duplicate_pairs"])
    rows.append([label, "exact", num(ex["pairs_compared"]),
                 num(ex["picard_duplicate_rate_pct"]),
                 num(ex["sequence_duplicate_rate_pct"]),
                 num(ex["sequence_duplicate_pairs"]),
                 num(ex["sequence_extra_copies_flagged_by_picard_pct"]),
                 round(found, 4)])
    for r in tsv(os.path.join(DATA, f"picard_groups_{k}.tsv")):
        rows.append([label, r["label"]] + [num(r[m]) for m, _ in GROUPS])
sheet(wb, "Mismatches allowed", ["Library", "Counting"] +
      [h for _, h in GROUPS], rows)
for k, label in SETS:
    mm = tsv(os.path.join(DATA, f"picard_set_mismatches_{k}.tsv"))
    sheet(wb, f"Picard sets, {k}",
          ["Mismatching bases (first 50 bases of both mates)",
           "Pairs, mate order", "Pairs (%), mate order",
           "Pairs, either mate order", "Pairs (%), either mate order"],
          [[r["mismatches"], num(r["pairs"]), num(r["pct"]),
            num(r["pairs_either_order"]), num(r["pct_either_order"])]
           for r in mm if r["mismatches"] not in ("length differs", "with_N",
                                                  "mates_swapped")])
wb.save(os.path.join(OUT, "S4_Table.xlsx"))

# S5 Table: comparisons with other programs
wb = Workbook()
opt = tsv(os.path.join(DATA, "optical_validation", "summary.tsv"))
sheet(wb, "Optical, all pairs",
      ["Run", "Distance (pixels)", "Pairs examined by Picard",
       "Picard duplicates", "Picard optical duplicates",
       "Mismatches allowed", "Pairs given to falcoAmadeus",
       "Pairs left out (would not align)", "falcoAmadeus duplicates",
       "falcoAmadeus optical duplicates", "falcoAmadeus wall time (s)",
       "falcoAmadeus peak memory (kB)"],
      [[r["dataset"], int(r["distance"]), int(r["picard_pairs_examined"]),
        int(r["picard_duplicates"]), int(r["picard_optical"]),
        int(r["mismatches"]), int(r["falco_pairs"]),
        int(r["falco_left_out"]), int(r["falco_duplicates"]),
        int(r["falco_optical"]), num(r["falco_seconds"]),
        int(r["falco_peak_kb"])] for r in opt], first=True)
OPT_KEYS = [("pairs_picard_examines", "Pairs examined by Picard"),
            ("falco_duplicates_on_these", "falcoAmadeus duplicates"),
            ("falco_optical_on_these_100",
             "falcoAmadeus optical duplicates, 100 pixels"),
            ("falco_optical_on_these_2500",
             "falcoAmadeus optical duplicates, 2500 pixels"),
            ("close_copies_same_picard_coordinates",
             "Close copies (100 pixels) in the same Picard set"),
            ("close_copies_other_picard_coordinates",
             "Close copies (100 pixels) at other Picard coordinates"),
            ("of_which_both_mapq_20_or_more",
             "  of which both pairs mapping quality >= 20")]
pp = {}
with open(os.path.join(DATA, "optical_validation",
                       "on_picard_pairs.tsv")) as f:
    for line in f:
        ds, m, key, val = line.rstrip("\n").split("\t")
        pp[(ds, m, key.strip())] = int(val)
picard = {(r["dataset"], r["distance"]): int(r["picard_optical"])
          for r in opt}
rows = []
for ds in sorted({k[0] for k in pp}):
    for m in sorted({k[1] for k in pp if k[0] == ds}):
        rows.append([ds, int(m[2:]), picard[(ds, "100")],
                     picard[(ds, "2500")]] +
                    [pp[(ds, m, k)] for k, _ in OPT_KEYS])
sheet(wb, "Optical, Picard pairs",
      ["Run", "Mismatches allowed", "Picard optical duplicates, 100 pixels",
       "Picard optical duplicates, 2500 pixels"] + [h for _, h in OPT_KEYS],
      rows)
trim = tsv(os.path.join(DATA, "trim_preview_validation",
                        "cutadapt_comparison.tsv"))
sheet(wb, "Trimming preview",
      ["File (first 200,000 reads)", "Cutoff", "Value (%)", "Preview",
       "cutadapt 5.2", "Equal"],
      [[r["file"], r["cutoff"], r["field"], num(r["preview"]),
        num(r["cutadapt"]), r["equal"]] for r in trim])
cx = tsv(os.path.join(DATA, "complexity_insert_validation", "summary.tsv"))
sheet(wb, "Complexity, insert size",
      ["Run", "Read pairs", "Comparison", "Values compared", "Equal"],
      [[r["dataset"], int(r["pairs"]), r["check"],
        int(r["values_compared"]), r["equal"]] for r in cx])
truth = {r["name"]: r for r in tsv(os.path.join(
    DATA, "screen_validation", "mixture_truth.tsv"))}
mix = []
with open(os.path.join(DATA, "screen_validation",
                       "mixture_screen.tsv")) as f:
    on = False
    for line in f:
        p = line.rstrip("\n").split("\t")
        if p[0] == "#Reference":
            on = True
            continue
        if on and len(p) >= 2:
            t = truth.get(p[0])
            mix.append([p[0], num(t["fraction"]) * 100 if t else 0,
                        int(t["reads"]) if t else 0, num(p[1]),
                        num(p[2]) if len(p) > 2 and p[2] else ""])
sheet(wb, "Screen, simulated",
      ["Reference", "Reads drawn (%)", "Reads drawn",
       "k-mers found (%)", "k-mers found only here (%)"], mix)
pub = tsv(os.path.join(DATA, "screen_validation", "public_runs.tsv"))
sheet(wb, "Screen, public runs",
      ["Run", "Description", "Reference", "k-mers found (%)",
       "k-mers found only here (%)"],
      [[r["run"], r["description"], r["reference"], num(r["found_pct"]),
        num(r["only_pct"])] for r in pub])
man = tsv(os.path.join(DATA, "screen_validation",
                       "references_manifest.tsv"))
sheet(wb, "Screen panel",
      ["Reference", "File", "Source", "Release", "Sequences", "Bases",
       "SHA-256"],
      [[r["name"], r["file"], r["source_url_or_accessions"], r["release"],
        int(r["sequences"]), int(r["bases"]), r["sha256"]] for r in man])
wb.save(os.path.join(OUT, "S5_Table.xlsx"))
