#!/usr/bin/env python3
"""Figure panels and table values for the falcoAmadeus manuscript.

Reads the measurements in paper/data and writes one PDF per figure panel
(no panel letters; the figures are assembled in LaTeX) plus TSV files with
the values quoted in the tables.

Usage: make_figures.py <data_dir> <out_dir>
"""
import csv
import os
import statistics
import sys
from collections import defaultdict

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.ticker import (FuncFormatter, LogLocator,  # noqa: E402
                               NullFormatter)

DATA, OUT = sys.argv[1], sys.argv[2]
os.makedirs(OUT, exist_ok=True)

plt.rcParams.update({
    "font.family": "Liberation Sans",  # metric-compatible with Arial (PLOS)
    "font.size": 8,
    "axes.titlesize": 8,
    "axes.labelsize": 8,
    "legend.fontsize": 7,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "pdf.fonttype": 42,
})

COL = {
    "fa": "#1f4e8c",       # falcoAmadeus
    "fa_pe": "#1f4e8c",
    "f2": "#d9822b",       # falco 2.0.3
    "f2o": "#d9822b",
    "f1": "#7f7f7f",       # falco 1.3.2
    "fastqc": "#3b3b3b",
}
LABEL = {
    "fa": "falcoAmadeus",
    "fa_pe": "falcoAmadeus paired mode",
    "fa_two": "falcoAmadeus, two files",
    "f2": "falco 2.0.3",
    "f2o": "falco 2.0.3 --orig-dups",
    "f1": "falco 1.3.2",
    "fastqc": "FastQC 0.12.1",
}
PANEL_W, PANEL_H = 2.6, 2.1


def read_tsv(name):
    with open(os.path.join(DATA, name)) as f:
        return list(csv.DictReader(f, delimiter="\t"))


def plain_log_y(ax):
    ax.yaxis.set_major_locator(LogLocator(base=10, subs=(1.0, 2.0, 5.0)))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
    ax.yaxis.set_minor_formatter(NullFormatter())


def legend_panel(ax_src, name, ncol):
    handles, labels = ax_src.get_legend_handles_labels()
    rows = -(-len(labels) // ncol)
    fig = plt.figure(figsize=(2 * PANEL_W, 0.2 + 0.17 * rows))
    fig.legend(handles, labels, loc="center", ncol=ncol, frameon=False)
    save(fig, name)


def save(fig, name):
    fig.savefig(os.path.join(OUT, name), bbox_inches="tight")
    plt.close(fig)


def medians(rows, key_fields):
    acc = defaultdict(lambda: ([], []))
    ok = defaultdict(lambda: True)
    for r in rows:
        k = tuple(r[f] for f in key_fields)
        acc[k][0].append(float(r["wall_s"]))
        acc[k][1].append(int(r["rss_kb"]) / 1024.0)
        ok[k] = ok[k] and r["exit"] == "0"
    return {k: (statistics.median(v[0]), max(v[1]), ok[k])
            for k, v in acc.items()}


# ---------------------------------------------------------------------------
# Single-end benchmark (whole machine, warm page cache)
se = (read_tsv("se_falcoamadeus.tsv") + read_tsv("se_falco2.tsv"))
parent_path = os.path.join(DATA, "se_parent_tools.tsv")
if os.path.exists(parent_path):
    se += read_tsv("se_parent_tools.tsv")
SE = medians(se, ("file", "shape", "cache", "threads", "tool"))
THREADS_WIDE = [1, 4, 8, 16, 32]
SE_FILES = [
    ("SRR14777166.fastq", "fig1a.pdf", "FASTQ"),
    ("SRR14777166.fastq.gz", "fig1b.pdf", "gzip FASTQ"),
    ("SRR14777166.bgzf.fastq.gz", "fig1c.pdf", "BGZF FASTQ"),
    ("SRR14777166.sam", "fig1d.pdf", "SAM"),
    ("SRR14777166.bam", "fig1e.pdf", "BAM"),
]


def runtime_panel(table, file_key, name, title, tools, key_fn, ref_tools):
    fig, ax = plt.subplots(figsize=(PANEL_W, PANEL_H))
    for tool, ls in tools:
        xs, ys = [], []
        for t in THREADS_WIDE:
            v = table.get(key_fn(file_key, t, tool))
            if v and v[2]:
                xs.append(t)
                ys.append(v[0])
        if xs:
            ax.plot(xs, ys, ls, color=COL[tool.split(":")[0]], marker="o",
                    ms=3, lw=1.2, label=LABEL[tool.split(":")[-1]])
    for tool, mk in ref_tools:
        v = table.get(key_fn(file_key, 1, tool))
        if v and v[2]:
            ax.plot([1], [v[0]], mk, color=COL[tool], ms=5, ls="none",
                    label=LABEL[tool])
    ax.set_xscale("log", base=2)
    ax.set_yscale("log")
    ax.set_xticks(THREADS_WIDE)
    ax.set_xticklabels([str(t) for t in THREADS_WIDE])
    ax.set_xlabel("Threads")
    ax.set_ylabel("Wall time (s)")
    ax.set_title(title)
    return fig, ax


for fkey, pname, title in SE_FILES:
    fig, ax = runtime_panel(
        SE, fkey, pname, title,
        [("fa", "-"), ("f2", "-"), ("f2o", "--")],
        lambda f, t, tool: (f, "wide", "warm", str(t), tool),
        [("f1", "s"), ("fastqc", "^")])
    for line in ax.get_lines():
        if line.get_label() == LABEL["f2o"]:
            line.set_marker("o")
            line.set_markerfacecolor("white")
    plain_log_y(ax)
    if pname == "fig1a.pdf":
        legend_panel(ax, "fig1legend.pdf", 3)
    save(fig, pname)

# memory and cold-cache values for the tables
with open(os.path.join(OUT, "table_se.tsv"), "w") as out:
    out.write("file\tshape\tcache\tthreads\ttool\twall_s\trss_mb\tok\n")
    for k, v in sorted(SE.items()):
        out.write("\t".join(k) + f"\t{v[0]:.2f}\t{v[1]:.0f}\t{int(v[2])}\n")

# ---------------------------------------------------------------------------
# Paired-end benchmark (10M HG002 pairs)
pe = (read_tsv("pe_competitors_and_twofile.tsv")
      + read_tsv("pe_falcoamadeus_paired.tsv"))
PE = medians(pe, ("format", "shape", "threads", "tool"))
for fmt, pname, title in [("fastq", "fig2a.pdf", "FASTQ"),
                          ("fastq.gz", "fig2b.pdf", "gzip FASTQ")]:
    fig, ax = plt.subplots(figsize=(PANEL_W, PANEL_H))
    series = [("fa_pe", "-", "fa_pe"), ("fa", ":", "fa_two"),
              ("f2", "-", "f2"), ("f2o", "--", "f2o")]
    for tool, ls, lab in series:
        xs, ys = [], []
        for t in THREADS_WIDE:
            v = PE.get((fmt, "wide", str(t), tool))
            if v and v[2]:
                xs.append(t)
                ys.append(v[0])
        mfc = "white" if tool in ("f2o", "fa") else None
        ax.plot(xs, ys, ls, color=COL[tool], marker="o", ms=3, lw=1.2,
                mfc=mfc, label=LABEL[lab])
    for tool, mk in [("f1", "s"), ("fastqc", "^")]:
        v = PE.get((fmt, "wide", "1", tool))
        if v and v[2]:
            ax.plot([1], [v[0]], mk, color=COL[tool], ms=5, ls="none",
                    label=LABEL[tool])
    ax.set_xscale("log", base=2)
    ax.set_yscale("log")
    ax.set_xticks(THREADS_WIDE)
    ax.set_xticklabels([str(t) for t in THREADS_WIDE])
    ax.set_xlabel("Threads")
    ax.set_ylabel("Wall time (s)")
    ax.set_title(title)
    plain_log_y(ax)
    if pname == "fig2a.pdf":
        legend_panel(ax, "fig2legend.pdf", 3)
    save(fig, pname)

with open(os.path.join(OUT, "table_pe.tsv"), "w") as out:
    out.write("format\tshape\tthreads\ttool\twall_s\trss_mb\tok\n")
    for k, v in sorted(PE.items()):
        out.write("\t".join(k) + f"\t{v[0]:.2f}\t{v[1]:.0f}\t{int(v[2])}\n")

# ---------------------------------------------------------------------------
# Paired-end duplication on simulated pairs
sim = read_tsv("duplication_simulation.tsv")
both = [r for r in sim if r["scenario"] == "both_mates"]
fig, ax = plt.subplots(figsize=(PANEL_W, PANEL_H))
ax.plot([0, 65], [0, 65], color="#bbbbbb", lw=0.8, zorder=0)
place_style = {"late": ("v", "--"), "early": ("^", ":"), "shuffled": ("o", "-.")}
for pl, (mk, ls) in place_style.items():
    rows = sorted((r for r in both if r["placement"] == pl),
                  key=lambda r: float(r["true_dup"]))
    ax.plot([100 * float(r["true_dup"]) for r in rows],
            [float(r["paired_fastqc"]) for r in rows], ls, marker=mk, ms=3,
            lw=1, color=COL["f2"], mfc="white",
            label=f"FastQC-style, copies {pl}")
rows = sorted(both, key=lambda r: float(r["true_dup"]))
ax.plot([100 * float(r["true_dup"]) for r in rows],
        [float(r["paired_exact"]) for r in rows], "o", ms=3.5,
        color=COL["fa"], label="Exact, all placements")
ax.set_xlabel("True duplicate pairs (%)")
ax.set_ylabel("Reported duplicate pairs (%)")
ax.set_xlim(-2, 65)
ax.set_ylim(-2, 70)
legend_panel(ax, "fig3legend.pdf", 2)
save(fig, "fig3a.pdf")

shared = sorted((r for r in sim if r["scenario"] == "shared_r1"),
                key=lambda r: float(r["true_dup"]))
fig, ax = plt.subplots(figsize=(PANEL_W, PANEL_H))
x = [100 * float(r["true_dup"]) for r in shared]
ax.plot(x, [float(r["r1_exact"]) for r in shared], "s-", ms=3, lw=1.2,
        color=COL["f1"], label="Read 1 alone")
ax.plot(x, [float(r["paired_exact"]) for r in shared], "o-", ms=3, lw=1.2,
        color=COL["fa"], label="Read pairs")
ax.set_xlabel("Pairs repeating read 1 only (%)")
ax.set_ylabel("Reported duplicates (%)")
ax.set_xlim(-2, 65)
ax.set_ylim(-2, 65)
ax.legend(frameon=False, loc="upper left")
save(fig, "fig3b.pdf")

# ---------------------------------------------------------------------------
# Assay-profile validation on 120 public runs
grades = read_tsv("assay_validation_grades.tsv")
SKIP = {"Basic Statistics", "Per tile sequence quality"}
TARGET = {
    "adapter": "Adapter Content",
    "qualtail": "Per base sequence quality",
    "midbias": "Per base sequence content",
    "gccontam": "Per sequence GC content",
    "dupspike": "Sequence Duplication Levels",
    "ncontent": "Per base N content",
}
flag = defaultdict(lambda: [0, 0])     # (assay, kind) -> flagged, total
det = defaultdict(lambda: [0, 0])      # (defect, kind) -> detected, total
mod_flag = defaultdict(lambda: [0, 0, 0])  # (assay, module) -> gen, prof, n
for r in grades:
    if r["module"] == "RUN_FAILED":
        continue
    kind = "generic" if r["profile"] == "generic" else "profile"
    flagged = r["grade"] in ("WARN", "FAIL")
    if r["condition"] == "clean" and r["module"] not in SKIP:
        f = flag[(r["assay"], kind)]
        f[0] += flagged
        f[1] += 1
        m = mod_flag[(r["assay"], r["module"])]
        m[0 if kind == "generic" else 1] += flagged
        if kind == "generic":
            m[2] += 1
    elif TARGET.get(r["condition"]) == r["module"]:
        d = det[(r["condition"], kind)]
        d[0] += flagged
        d[1] += 1

ASSAY_ORDER = ["wgs", "exome", "chipseq", "atac", "rnaseq", "smallrna",
               "amplicon16s", "wgbs", "emseq", "rrbs"]
ASSAY_LABEL = {"wgs": "WGS", "exome": "Exome", "chipseq": "ChIP-seq",
               "atac": "ATAC-seq", "rnaseq": "RNA-seq",
               "smallrna": "Small RNA", "amplicon16s": "16S",
               "wgbs": "WGBS", "emseq": "EM-seq", "rrbs": "RRBS"}
fig, ax = plt.subplots(figsize=(2 * PANEL_W, PANEL_H))
w = 0.38
xs = range(len(ASSAY_ORDER))
gen = [100 * flag[(a, "generic")][0] / flag[(a, "generic")][1]
       for a in ASSAY_ORDER]
pro = [100 * flag[(a, "profile")][0] / flag[(a, "profile")][1]
       for a in ASSAY_ORDER]
ax.bar([x - w / 2 for x in xs], gen, w, color="#b0b0b0", label="Generic grading")
ax.bar([x + w / 2 for x in xs], pro, w, color=COL["fa"], label="Assay profile")
ax.set_xticks(list(xs))
ax.set_xticklabels([ASSAY_LABEL[a] for a in ASSAY_ORDER], rotation=40,
                   ha="right")
ax.set_ylabel("Module checks flagged (%)")
ax.legend(frameon=False, loc="upper left", ncol=2)
save(fig, "fig4a.pdf")

DEFECT_ORDER = ["adapter", "qualtail", "gccontam", "dupspike", "ncontent",
                "midbias"]
DEFECT_LABEL = {"adapter": "Adapter", "qualtail": "Quality tail",
                "gccontam": "GC contamination", "dupspike": "Duplicate spike",
                "midbias": "Composition bias", "ncontent": "N block"}
fig, ax = plt.subplots(figsize=(PANEL_W, PANEL_H))
xs = range(len(DEFECT_ORDER))
gen = [100 * det[(d, "generic")][0] / det[(d, "generic")][1]
       for d in DEFECT_ORDER]
pro = [100 * det[(d, "profile")][0] / det[(d, "profile")][1]
       for d in DEFECT_ORDER]
ax.bar([x - w / 2 for x in xs], gen, w, color="#b0b0b0", label="Generic grading")
ax.bar([x + w / 2 for x in xs], pro, w, color=COL["fa"], label="Assay profile")
ax.set_xticks(list(xs))
ax.set_xticklabels([DEFECT_LABEL[d] for d in DEFECT_ORDER], rotation=40,
                   ha="right")
ax.set_ylabel("Runs with defect detected (%)")
ax.set_ylim(0, 115)
ax.set_yticks([0, 25, 50, 75, 100])
save(fig, "fig4b.pdf")

with open(os.path.join(OUT, "table_assay_modules.tsv"), "w") as out:
    out.write("assay\tmodule\tgeneric_flagged\tprofile_flagged\truns\n")
    for (a, m), v in sorted(mod_flag.items()):
        out.write(f"{a}\t{m}\t{v[0]}\t{v[1]}\t{v[2]}\n")
with open(os.path.join(OUT, "table_assay_detection.tsv"), "w") as out:
    out.write("defect\tgeneric_detected\tprofile_detected\truns\n")
    for d in DEFECT_ORDER:
        out.write(f"{d}\t{det[(d, 'generic')][0]}\t{det[(d, 'profile')][0]}"
                  f"\t{det[(d, 'generic')][1]}\n")
