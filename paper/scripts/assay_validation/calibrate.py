#!/usr/bin/env python3
"""Calibrates the assay-profile limits on runs that are NOT the validation
runs (an independent calibration set from other BioProjects).

Usage: calibrate.py <runs.tsv> <fastq_dir> <falcoAmadeus> <out_dir> [jobs]

For every clean calibration run: the generic grading metrics
(FA_GRADE_METRICS: per-base maximum difference, GC deviation, deduplicated
percentage, adapter maximum), the per-position composition (--nogroup), and
the conversion estimate of the WGBS profile (on every run, so that
non-bisulfite runs show where unconverted libraries lie). Runs whose read 1 is shorter than 20 bases (a cell barcode, not genomic
sequence) are excluded and listed in excluded.tsv. Writes metrics.tsv
and limits.tsv:
- atac_mask: the largest position (from the 5' end, within 25) whose
  composition differs by more than 10 % in any ATAC run;
- gc_fail per relaxed assay: 1.25 x its largest clean GC deviation (rounded
  up to 5), at least the generic FAIL limit 30;
- dup_fail per relaxed assay: 0.8 x its lowest clean deduplicated percentage
  (rounded down to 1), at most the generic FAIL limit 50, and 0 when that
  percentage is below 2 % (no usable limit);
- conversion: the largest estimate of the bisulfite/EM runs and the
  smallest of the other runs.
"""
import concurrent.futures as cf
import csv
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile

RELAX_GC = ["exome", "amplicon16s", "rnaseq", "smallrna", "wgbs", "emseq",
            "rrbs"]
RELAX_DUP = ["exome", "amplicon16s", "rnaseq", "smallrna", "atac", "chipseq",
             "rrbs"]
METHYL = ["wgbs", "emseq", "rrbs"]


def run_fa(fa, fq, od, extra, env=None):
    p = subprocess.run([fa, "-q", "-t", "2", "-o", od, *extra, fq],
                       check=True, stdout=subprocess.DEVNULL,
                       stderr=subprocess.PIPE, text=True, env=env)
    return p.stderr


def one(assay, run, fqdir, fa, out):
    fq = os.path.join(fqdir, run + ".fastq")
    tmp = tempfile.mkdtemp(dir=out)
    try:
        env = dict(os.environ, FA_GRADE_METRICS="1")
        err = run_fa(fa, fq, os.path.join(tmp, "g"), [], env)
        met = {}
        for line in err.splitlines():
            if line.startswith("FA_METRIC\t"):
                _, m, v = line.split("\t")
                met[m] = float(v)
        run_fa(fa, fq, os.path.join(tmp, "n"), ["--nogroup"])
        rows, inside = [], False
        for line in open(os.path.join(tmp, "n", "fastqc_data.txt")):
            if line.startswith(">>Per base sequence content"):
                inside = True
            elif inside and line.startswith(">>END_MODULE"):
                break
            elif inside and not line.startswith("#"):
                x = line.rstrip("\n").split("\t")
                rows.append([float(v) for v in x[1:5]])  # G A T C
        dev = [max(abs(r[1] - r[2]), abs(r[0] - r[3])) for r in rows]
        bias5 = max([i + 1 for i in range(min(25, len(dev)))
                     if dev[i] > 10.0], default=0)
        run_fa(fa, fq, os.path.join(tmp, "w"), ["--assay", "wgbs"])
        conv = float("nan")
        for line in open(os.path.join(tmp, "w", "fastqc_data.txt")):
            m = re.search(r"conversion estimate: unconverted C outside "
                          r"CpG ([0-9.]+)%", line)
            if m:
                conv = float(m.group(1))
        return (assay, run, met.get("Per base sequence content", math.nan),
                met.get("Per sequence GC content", math.nan),
                met.get("Sequence Duplication Levels", math.nan),
                met.get("Adapter Content", math.nan), bias5, conv)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def main():
    runs_tsv, fqdir, fa, out = sys.argv[1:5]
    jobs = int(sys.argv[5]) if len(sys.argv) > 5 else 8
    os.makedirs(out, exist_ok=True)
    runs, excluded = [], []
    for r in csv.DictReader(open(runs_tsv), delimiter="\t"):
        fq = os.path.join(fqdir, r["run"] + ".fastq")
        if not os.path.exists(fq):
            continue
        with open(fq) as f:
            f.readline()
            n = len(f.readline().strip())
        if n < 20:  # read 1 is a barcode, not genomic sequence
            excluded.append((r["assay"], r["run"], f"read 1 is {n} bases"))
            continue
        runs.append((r["assay"], r["run"]))
    with open(os.path.join(out, "excluded.tsv"), "w") as f:
        f.write("assay\trun\treason\n")
        for x in excluded:
            f.write("\t".join(x) + "\n")
    with cf.ProcessPoolExecutor(jobs) as ex:
        res = list(ex.map(one, *zip(*[(a, r, fqdir, fa, out)
                                      for a, r in runs])))
    with open(os.path.join(out, "metrics.tsv"), "w") as f:
        f.write("assay\trun\tper_base_max_diff\tgc_deviation\tdedup_pct"
                "\tadapter_max\tbias_end_5p\tconversion_pct\n")
        for r in sorted(res):
            f.write("\t".join(str(x) for x in r) + "\n")
    by = {}
    for r in res:
        by.setdefault(r[0], []).append(r)
    lim = []
    lim.append(("atac_mask", "atac",
                max(r[6] for r in by.get("atac", [])) if "atac" in by
                else ""))
    for a in RELAX_GC:
        if a == "amplicon16s":
            # a few discrete products: the normal-curve deviation has no
            # meaning, no limit (as no per-base composition grade)
            lim.append(("gc_fail", a, 0))
            continue
        if a in by:
            g = max(r[3] for r in by[a])
            # never stricter than the generic FAIL limit (30)
            lim.append(("gc_fail", a, max(30, 5 * math.ceil(1.25 * g / 5))))
    for a in RELAX_DUP:
        if a in by:
            d = min(r[4] for r in by[a])
            # never stricter than the generic FAIL limit (50 % left)
            lim.append(("dup_fail", a,
                        min(50, math.floor(0.8 * d)) if d >= 2 else 0))
    meth = [r[7] for r in res if r[0] in METHYL]
    other = [r[7] for r in res if r[0] not in METHYL]
    lim.append(("conversion_max_bisulfite", "", max(meth) if meth else ""))
    lim.append(("conversion_min_other", "", min(other) if other else ""))
    with open(os.path.join(out, "limits.tsv"), "w") as f:
        f.write("limit\tassay\tvalue\n")
        for x in lim:
            f.write("\t".join(str(v) for v in x) + "\n")
    print(open(os.path.join(out, "limits.tsv")).read())


if __name__ == "__main__":
    main()
