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
    ("SRR14777166.fastq", "fig4a.pdf", "FASTQ"),
    ("SRR14777166.fastq.gz", "fig4b.pdf", "gzip FASTQ"),
    ("SRR14777166.bgzf.fastq.gz", "fig4c.pdf", "BGZF FASTQ"),
    ("SRR14777166.sam", "fig4d.pdf", "SAM"),
    ("SRR14777166.bam", "fig4e.pdf", "BAM"),
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
    if pname == "fig4a.pdf":
        legend_panel(ax, "fig4legend.pdf", 3)
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
for fmt, pname, title in [("fastq", "fig5a.pdf", "FASTQ"),
                          ("fastq.gz", "fig5b.pdf", "gzip FASTQ")]:
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
    if pname == "fig5a.pdf":
        legend_panel(ax, "fig5legend.pdf", 3)
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
legend_panel(ax, "fig1clegend.pdf", 2)
save(fig, "fig1c.pdf")

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
save(fig, "fig1d.pdf")

# Duplicate rates of alignment-free tools on a simulated library with known
# PCR duplication (simulate_duplicates.py, run_duplication_competitors.sh)
comp = read_tsv("duplication_competitors.tsv")
truth = {r["metric"]: float(r["value"])
         for r in read_tsv("duplication_competitors_truth.tsv")}


def comp_rate(tool, mode, setting):
    for r in comp:
        if (r["tool"], r["mode"], r["setting"]) == (tool, mode, setting):
            return float(r["duplicate_rate_pct"])
    raise KeyError((tool, mode, setting))


FA_DEF = "default (FastQC estimate; PE: exact pairs)"
FA_UMI = "2 mismatches, UMIs 1:8 (1 mismatch)"
COMP_GROUPS = [
    ("Read 1 alone", truth["se_pcr_duplicate_rate_pct"], [
        ("FastQC 0.12.1", "FastQC", "SE read 1", "default"),
        ("falco 1.3.2", "falco 1.3.2", "SE read 1", "default"),
        ("falco 2.0.3", "falco 2", "SE read 1", "default"),
        ("fastp", "fastp", "SE", "accuracy 3"),
        ("clumpify, 2 substitutions", "clumpify", "SE", "subs=2"),
        ("falcoAmadeus, 2 mismatches", "falcoAmadeus", "SE", "2 mismatches"),
        ("falcoAmadeus, 2 mismatches, UMI", "falcoAmadeus", "SE", FA_UMI),
    ]),
    ("Read pairs", truth["pe_pcr_duplicate_rate_pct"], [
        ("fastp", "fastp", "PE", "accuracy 3"),
        ("clumpify, 2 substitutions", "clumpify", "PE", "subs=2"),
        ("Picard EstimateLibraryComplexity", "Picard ELC", "PE",
         "MAX_DIFF_RATE 0.03"),
        ("Picard EstimateLibraryComplexity, UMI", "Picard ELC", "PE",
         "MAX_DIFF_RATE 0.03, UMIs as read barcode tags"),
        ("falcoAmadeus, exact", "falcoAmadeus", "PE", FA_DEF),
        ("falcoAmadeus, 2 mismatches", "falcoAmadeus", "PE", "2 mismatches"),
        ("falcoAmadeus, 2 mismatches, UMI", "falcoAmadeus", "PE", FA_UMI),
    ]),
]
fig, ax = plt.subplots(figsize=(2 * PANEL_W, 2.9))
y = 0
ticks, labels = [], []
for gname, true_rate, items in COMP_GROUPS:
    top = y
    for lab, tool, mode, setting in items:
        v = comp_rate(tool, mode, setting)
        ax.barh(y, v, 0.7, color=COL["fa"] if tool == "falcoAmadeus"
                else "#a8a8a8")
        ax.text(v + 0.4, y, f"{v:.1f}", va="center", fontsize=6.5)
        ticks.append(y)
        labels.append(lab)
        y += 1
    ax.plot([true_rate, true_rate], [top - 0.5, y - 0.5], "--", color="black",
            lw=0.9)
    ax.text(true_rate + 0.3, top - 0.75, f"{gname}: true {true_rate:.1f} %",
            fontsize=6.5, va="bottom")
    y += 1.2
ax.set_yticks(ticks)
ax.set_yticklabels(labels)
ax.invert_yaxis()
ax.set_xlabel("Reported duplicate rate (%)")
ax.set_xlim(0, 36)
save(fig, "fig1a.pdf")

# Error-aware groups of copies against Picard MarkDuplicates on the same
# pairs: Picard duplicates found and counted copies that Picard marks
DATASETS = [("drr041223", "RNA-seq (DRR041223)"),
            ("hg002", "WGS (HG002)")]
fig, ax = plt.subplots(figsize=(PANEL_W, PANEL_H))
mk = {"drr041223": "o", "hg002": "s"}
for key, lab in DATASETS:
    ex = {r["metric"]: float(r["value"])
          for r in read_tsv(f"picard_comparison_{key}.tsv")}
    found = [100 * ex["sequence_extra_copies_flagged_by_picard"]
             / ex["picard_duplicate_pairs"]]
    prec = [ex["sequence_extra_copies_flagged_by_picard_pct"]]
    for r in read_tsv(f"picard_groups_{key}.tsv"):
        found.append(float(r["picard_duplicates_found_pct"]))
        prec.append(float(r["extra_copies_flagged_by_picard_pct"]))
    ax.plot(range(3), found, "-", marker=mk[key], ms=3.5, lw=1.2,
            color=COL["fa"], label=f"{lab}, Picard duplicates found")
    ax.plot(range(3), prec, "--", marker=mk[key], ms=3.5, lw=1.2,
            color=COL["f2"], mfc="white",
            label=f"{lab}, counted copies marked by Picard")
ax.set_xticks(range(3))
ax.set_xticklabels(["0", "1", "2"])
ax.set_xlabel("Mismatches allowed")
ax.set_ylabel("Pairs (%)")
ax.set_ylim(0, 105)
legend_panel(ax, "fig1blegend.pdf", 1)
save(fig, "fig1b.pdf")

# ---------------------------------------------------------------------------
# Optical duplicates, library complexity and contamination screen
opt_sum = read_tsv("optical_validation/summary.tsv")
picard_opt = {(r["dataset"], r["distance"]): int(r["picard_optical"])
              for r in opt_sum}
on_pp = defaultdict(dict)
with open(os.path.join(DATA, "optical_validation/on_picard_pairs.tsv")) as f:
    for line in f:
        ds, m, key, val = line.rstrip("\n").split("\t")
        on_pp[(ds, m)][key.strip()] = int(val)
fig, ax = plt.subplots(figsize=(PANEL_W, PANEL_H))
OPT_SETS = [("hiseq2500_HG002", "HiSeq 2500 (HG002)", "o"),
            ("novaseq_ERR5654681", "NovaSeq 6000 (ERR5654681)", "s")]
for ds, lab, m in OPT_SETS:
    for dist, ls, mfc in (("100", "-", None), ("2500", "--", "white")):
        ys = [100 * on_pp[(ds, f"M={k}")][f"falco_optical_on_these_{dist}"]
              / picard_opt[(ds, dist)] for k in range(3)]
        ax.plot(range(3), ys, ls, marker=m, ms=3.5, lw=1.2, mfc=mfc,
                color=COL["fa"] if ds.startswith("hiseq") else COL["f2"],
                label=f"{lab}, {dist} pixels")
ax.axhline(100, color="#bbbbbb", lw=0.8, zorder=0)
ax.set_xticks(range(3))
ax.set_xticklabels(["0", "1", "2"])
ax.set_xlabel("Mismatches allowed")
ax.set_ylabel("Optical duplicates (% of Picard)")
ax.set_ylim(0, 160)
legend_panel(ax, "fig2alegend.pdf", 1)
save(fig, "fig2a.pdf")

fig, ax = plt.subplots(figsize=(PANEL_W, PANEL_H))
CX = [("DRR041223", "RNA-seq (DRR041223)", COL["fa"]),
      ("SRR12898355_10M", "WGS (HG002)", COL["f2"])]
for name, lab, c in CX:
    with open(os.path.join(DATA, "complexity_insert_validation",
                           f"{name}.yield.tsv")) as f:
        seq = [int(line.split("\t")[1]) for line in f
               if line.startswith("1x\t")][0]
    xs, ys = [], []
    with open(os.path.join(DATA, "complexity_insert_validation",
                           f"{name}.curve.tsv")) as f:
        next(f)
        for line in f:
            a, b = (float(v) for v in line.split("\t")[:2])
            if 0 < a <= 100 * seq:
                xs.append(a / 1e6)
                ys.append(b / 1e6)
    ax.plot([seq / 1e6], [ys[min(range(len(xs)),
                                 key=lambda i: abs(xs[i] - seq / 1e6))]],
            "o", color=c, ms=4)
    ax.plot(xs, ys, "-", color=c, lw=1.2, label=lab)
    ax.plot([xs[0], xs[-1]], [xs[0], xs[-1]], color="#bbbbbb", lw=0.8,
            zorder=0)
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel("Read pairs sequenced (millions)")
ax.set_ylabel("Distinct pairs (millions)")
ax.legend(frameon=False, loc="upper left")
save(fig, "fig2b.pdf")

mix_truth = {r["name"]: float(r["fraction"]) * 100
             for r in read_tsv("screen_validation/mixture_truth.tsv")}
mix = {}
with open(os.path.join(DATA, "screen_validation/mixture_screen.tsv")) as f:
    rows = [line.rstrip("\n").split("\t") for line in f]
start = [i for i, r in enumerate(rows) if r[0] == "#Reference"][0]
for r in rows[start + 1:]:
    if len(r) >= 3 and r[0] not in ("Several references", "No reference"):
        mix[r[0]] = float(r[1])
SCREEN_LABEL = {"human": "Human", "mouse": "Mouse", "rat": "Rat",
                "fly": "Fly", "worm": "Worm", "yeast": "Yeast",
                "arabidopsis": "Arabidopsis", "ecoli": "E. coli",
                "phix": "PhiX", "lambda": "Lambda",
                "human_mt": "Human mt", "mouse_mt": "Mouse mt",
                "human_rrna": "Human rRNA", "mycoplasma": "Mycoplasma",
                "vectors": "UniVec"}
order = sorted(mix, key=lambda k: (-mix_truth.get(k, 0), -mix[k]))
fig, ax = plt.subplots(figsize=(2 * PANEL_W, PANEL_H))
xs = range(len(order))
floor = 1e-3
ax.bar([x - 0.2 for x in xs], [max(mix_truth.get(k, 0), floor) for k in order],
       0.4, color="#a8a8a8", label="Simulated reads")
ax.bar([x + 0.2 for x in xs], [max(mix[k], floor) for k in order], 0.4,
       color=COL["fa"], label="Screen: k-mers found")
ax.set_yscale("log")
ax.set_ylim(floor, 200)
ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
ax.set_xticks(list(xs))
ax.set_xticklabels([SCREEN_LABEL[k] for k in order], rotation=40, ha="right")
ax.set_ylabel("Share (%)")
ax.legend(frameon=False, loc="upper right", ncol=2)
save(fig, "fig2c.pdf")

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
save(fig, "fig3a.pdf")

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
save(fig, "fig3b.pdf")

with open(os.path.join(OUT, "table_assay_modules.tsv"), "w") as out:
    out.write("assay\tmodule\tgeneric_flagged\tprofile_flagged\truns\n")
    for (a, m), v in sorted(mod_flag.items()):
        out.write(f"{a}\t{m}\t{v[0]}\t{v[1]}\t{v[2]}\n")
with open(os.path.join(OUT, "table_assay_detection.tsv"), "w") as out:
    out.write("defect\tgeneric_detected\tprofile_detected\truns\n")
    for d in DEFECT_ORDER:
        out.write(f"{d}\t{det[(d, 'generic')][0]}\t{det[(d, 'profile')][0]}"
                  f"\t{det[(d, 'generic')][1]}\n")
