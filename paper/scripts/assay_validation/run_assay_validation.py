#!/usr/bin/env python3
"""Grades of every assay profile on public runs of each assay type, clean
and with injected defects.

Usage: run_assay_validation.py <runs.tsv> <fastq_dir> <falcoAmadeus>
                               <inject_defect> <out_dir> [jobs]

runs.tsv: tab-separated with columns `assay` and `run`; <fastq_dir>/<run>.fastq
holds the reads. For every run:
- clean reads, graded with every profile (cross matrix), and once more with
  the generic profile and --nogroup, whose per-position base composition is
  kept for checking the profile masks;
- each defect of inject_defect, graded with the generic profile and with the
  run's own profile.
Writes grades.tsv (assay, run, condition, profile, module, grade),
annotations.tsv (the #Assay lines of fastqc_data.txt) and
per_base/<run>.tsv (position, G, A, T, C in percent).
"""
import concurrent.futures as cf
import csv
import os
import shutil
import subprocess
import sys
import tempfile

PROFILES = ["generic", "wgs", "exome", "amplicon", "16s", "rnaseq",
            "smallrna", "atac", "chipseq", "wgbs", "emseq", "rrbs"]
OWN = {"amplicon16s": "16s"}  # assay label -> profile, where they differ
DEFECTS = ["adapter", "qualtail", "midbias", "gccontam", "ncontent",
           "dupspike"]


def falco(fa, fq, outdir, profile, extra=()):
    cmd = [fa, "-q", "-t", "2", "-o", outdir, *extra]
    if profile != "generic":
        cmd += ["--assay", profile]
    subprocess.run(cmd + [fq], check=True, stdout=subprocess.DEVNULL,
                   stderr=subprocess.DEVNULL)
    grades = []
    with open(os.path.join(outdir, "summary.txt")) as f:
        for line in f:
            g, m = line.rstrip("\n").split("\t")[:2]
            grades.append((m, g))
    notes = []
    with open(os.path.join(outdir, "fastqc_data.txt")) as f:
        for line in f:
            if line.startswith("#Assay"):
                notes.append(line.rstrip("\n"))
    return grades, notes


def per_base(path):
    rows, inside = [], False
    with open(path) as f:
        for line in f:
            if line.startswith(">>Per base sequence content"):
                inside = True
                continue
            if inside:
                if line.startswith(">>END_MODULE"):
                    break
                if line.startswith("#"):
                    continue
                rows.append(line.rstrip("\n").split("\t"))
    return rows


def one_run(assay, run, fqdir, fa, inj, out):
    own = OWN.get(assay, assay)
    fq = os.path.join(fqdir, run + ".fastq")
    g_rows, n_rows = [], []
    tmp = tempfile.mkdtemp(dir=out)
    try:
        for p in PROFILES:
            od = os.path.join(tmp, "c_" + p)
            grades, notes = falco(fa, fq, od, p)
            g_rows += [(assay, run, "clean", p, m, g) for m, g in grades]
            n_rows += [(assay, run, "clean", p, x) for x in notes]
        od = os.path.join(tmp, "nogroup")
        falco(fa, fq, od, "generic", ["--nogroup"])
        with open(os.path.join(out, "per_base", run + ".tsv"), "w") as f:
            f.write("position\tG\tA\tT\tC\n")
            for r in per_base(os.path.join(od, "fastqc_data.txt")):
                f.write("\t".join(r[:5]) + "\n")
        for i, d in enumerate(DEFECTS):
            dfq = os.path.join(tmp, d + ".fastq")
            seed = str(sum(map(ord, run + d)) * 7919 + i)  # deterministic
            subprocess.run([inj, d, seed, fq, dfq], check=True)
            for p in dict.fromkeys(["generic", own]):
                grades, notes = falco(fa, dfq, os.path.join(tmp, d + p), p)
                g_rows += [(assay, run, d, p, m, g) for m, g in grades]
                n_rows += [(assay, run, d, p, x) for x in notes]
            os.remove(dfq)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    return g_rows, n_rows


def main():
    runs_tsv, fqdir, fa, inj, out = sys.argv[1:6]
    jobs = int(sys.argv[6]) if len(sys.argv) > 6 else 8
    os.makedirs(os.path.join(out, "per_base"), exist_ok=True)
    runs = [(r["assay"], r["run"])
            for r in csv.DictReader(open(runs_tsv), delimiter="\t")]
    with open(os.path.join(out, "grades.tsv"), "w") as fg, \
         open(os.path.join(out, "annotations.tsv"), "w") as fn, \
         cf.ProcessPoolExecutor(jobs) as ex:
        fg.write("assay\trun\tcondition\tprofile\tmodule\tgrade\n")
        fn.write("assay\trun\tcondition\tprofile\tannotation\n")
        futs = {ex.submit(one_run, a, r, fqdir, fa, inj, out): (a, r)
                for a, r in runs}
        for k, fu in enumerate(cf.as_completed(futs), 1):
            g_rows, n_rows = fu.result()
            fg.writelines("\t".join(x) + "\n" for x in g_rows)
            fn.writelines("\t".join(x) + "\n" for x in n_rows)
            print(f"{k}/{len(runs)} {futs[fu][0]} {futs[fu][1]}",
                  file=sys.stderr, flush=True)


if __name__ == "__main__":
    main()
