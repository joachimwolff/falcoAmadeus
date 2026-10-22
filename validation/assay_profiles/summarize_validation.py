#!/usr/bin/env python3
"""Summarize run_validation.sh output.

For each assay and profile (generic vs. the assay's own profile):

  clean     fraction of runs with each module flagged (WARN or FAIL) on the
            unmodified reads: the false-alarm rate on public libraries
  defects   fraction of runs in which the module targeted by an injected
            defect is flagged: the detection rate

Usage: summarize_validation.py validation.tsv [out_prefix]
Writes <out_prefix>_clean.tsv and <out_prefix>_detection.tsv and prints both.
"""
import csv
import sys
from collections import defaultdict

TARGET = {
    "adapter": "Adapter Content",
    "qualtail": "Per base sequence quality",
    "midbias": "Per base sequence content",
    "gccontam": "Per sequence GC content",
    "dupspike": "Sequence Duplication Levels",
    "ncontent": "Per base N content",
}


def main():
    path = sys.argv[1]
    prefix = sys.argv[2] if len(sys.argv) > 2 else path.rsplit(".", 1)[0]
    grades = {}
    runs = defaultdict(set)
    failed = []
    with open(path) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if row["module"] == "RUN_FAILED":
                failed.append(row)
                continue
            key = (row["assay"], row["run"], row["condition"], row["profile"],
                   row["module"])
            grades[key] = row["grade"]
            runs[row["assay"]].add(row["run"])

    def is_generic(profile):
        return profile == "generic"

    # clean: per assay, profile kind, module -> flagged count / runs graded
    clean = defaultdict(lambda: [0, 0])
    detect = defaultdict(lambda: [0, 0])
    for (assay, run, cond, prof, module), g in grades.items():
        kind = "generic" if is_generic(prof) else "profile"
        flagged = g in ("WARN", "FAIL")
        if cond == "clean":
            c = clean[(assay, kind, module)]
            c[0] += flagged
            c[1] += 1
        elif TARGET.get(cond) == module:
            d = detect[(assay, kind, cond)]
            d[0] += flagged
            d[1] += 1

    modules = sorted({k[2] for k in clean})
    with open(prefix + "_clean.tsv", "w") as out:
        out.write("assay\tmodule\tgeneric_flagged\tprofile_flagged\truns\n")
        for assay in sorted(runs):
            for m in modules:
                g = clean.get((assay, "generic", m), [0, 0])
                p = clean.get((assay, "profile", m), [0, 0])
                n = max(g[1], p[1])
                if n:
                    out.write(f"{assay}\t{m}\t{g[0]}\t{p[0]}\t{n}\n")
    with open(prefix + "_detection.tsv", "w") as out:
        out.write("assay\tdefect\ttarget_module\tgeneric_detected\t"
                  "profile_detected\truns\n")
        for assay in sorted(runs):
            for d, m in TARGET.items():
                g = detect.get((assay, "generic", d), [0, 0])
                p = detect.get((assay, "profile", d), [0, 0])
                n = max(g[1], p[1])
                if n:
                    out.write(f"{assay}\t{d}\t{m}\t{g[0]}\t{p[0]}\t{n}\n")
    for suffix in ("_clean.tsv", "_detection.tsv"):
        print(open(prefix + suffix).read())
    if failed:
        print(f"{len(failed)} grading runs failed", file=sys.stderr)


if __name__ == "__main__":
    main()
