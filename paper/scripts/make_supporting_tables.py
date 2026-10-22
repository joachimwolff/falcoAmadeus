#!/usr/bin/env python3
"""Supporting tables S1 to S3 of the falcoAmadeus manuscript as Excel files.

Usage: make_supporting_tables.py <data_dir> <panels_dir> <out_dir>
  panels_dir holds table_se.tsv, table_pe.tsv, table_assay_modules.tsv and
  table_assay_detection.tsv written by make_figures.py.
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
         "midbias": "Composition bias (Per base sequence content)"}
sheet(wb, "Injected defects",
      ["Defect (targeted module)", "Runs detected, generic",
       "Runs detected, profile", "Runs"],
      [[LABEL[r["defect"]], int(r["generic_detected"]),
        int(r["profile_detected"]), int(r["runs"])]
       for r in det if r["defect"] in LABEL])
wb.save(os.path.join(OUT, "S3_Table.xlsx"))
