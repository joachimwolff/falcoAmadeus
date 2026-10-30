#!/usr/bin/env python3
"""Tables of the assay-profile validation (run_assay_validation.py output).

Usage: analyze_assay_validation.py <out_dir>

Writes into <out_dir>:
- clean_by_module.tsv: per assay and module, runs flagged (WARN or FAIL)
  and failed with the generic profile and with the assay's own profile;
- defects.tsv: per assay and injected defect, runs whose target module is
  flagged, generic and own profile;
- defects_escalation.tsv: runs whose target module grades worse with the
  defect than the same run's clean reads under the same profile (a defect
  graded WARN on a library already WARN is not seen);
- cross_profiles.tsv: per assay of the data and profile, mean number of
  flagged and failed modules per run;
- masks.tsv: per run of a masked profile, how far the per-position
  composition bias reaches from the 5' (and 3') end, against the mask;
- conversion.tsv: the bisulfite/EM conversion annotations per assay and
  profile.
"""
import collections
import csv
import os
import re
import statistics as st
import sys

OUT = sys.argv[1]
OWN = {"amplicon16s": "16s"}
TARGET = {"adapter": "Adapter Content",
          "qualtail": "Per base sequence quality",
          "midbias": "Per base sequence content",
          "gccontam": "Per sequence GC content",
          "ncontent": "Per base N content",
          "dupspike": "Sequence Duplication Levels"}
# mask_per_base, mask_3p, methylation rule (FalcoConfig.cpp)
MASK = {"rnaseq": (13, 0, False), "atac": (11, 0, False),
        "rrbs": (3, 0, True), "emseq": (10, 10, True), "wgbs": (0, 0, True)}
WARN, FAIL = 10.0, 20.0  # Per base sequence content limits

grades = list(csv.DictReader(open(os.path.join(OUT, "grades.tsv")),
                             delimiter="\t"))
by = collections.defaultdict(dict)  # (assay, run, cond, profile) -> module
for g in grades:
    by[(g["assay"], g["run"], g["condition"], g["profile"])][g["module"]] = \
        g["grade"]
assays = sorted({g["assay"] for g in grades})
modules = list(dict.fromkeys(g["module"] for g in grades))
runs = collections.defaultdict(list)
for a, r, c, p in by:
    if c == "clean" and p == "generic":
        runs[a].append(r)


def flagged(gr):
    return gr in ("WARN", "FAIL")


# 1. clean runs, per module
with open(os.path.join(OUT, "clean_by_module.tsv"), "w") as f:
    f.write("assay\tprofile\tmodule\truns\tgeneric_flagged\tprofile_flagged"
            "\tgeneric_failed\tprofile_failed\n")
    for a in assays:
        own = OWN.get(a, a)
        for m in modules:
            gen = [by[(a, r, "clean", "generic")].get(m) for r in runs[a]]
            pro = [by[(a, r, "clean", own)].get(m) for r in runs[a]]
            f.write(f"{a}\t{own}\t{m}\t{len(runs[a])}\t"
                    f"{sum(map(flagged, gen))}\t{sum(map(flagged, pro))}\t"
                    f"{gen.count('FAIL')}\t{pro.count('FAIL')}\n")

# 2. injected defects
with open(os.path.join(OUT, "defects.tsv"), "w") as f:
    f.write("assay\tprofile\tdefect\ttarget_module\truns\tgeneric_flagged"
            "\tprofile_flagged\tgeneric_failed\tprofile_failed"
            "\tclean_profile_flagged\n")
    for a in assays:
        own = OWN.get(a, a)
        for d, m in TARGET.items():
            gen = [by[(a, r, d, "generic")].get(m) for r in runs[a]]
            pro = [by[(a, r, d, own)].get(m) for r in runs[a]]
            cln = [by[(a, r, "clean", own)].get(m) for r in runs[a]]
            f.write(f"{a}\t{own}\t{d}\t{m}\t{len(runs[a])}\t"
                    f"{sum(map(flagged, gen))}\t{sum(map(flagged, pro))}\t"
                    f"{gen.count('FAIL')}\t{pro.count('FAIL')}\t"
                    f"{sum(map(flagged, cln))}\n")

# 2b. a defect is seen only if it raises the grade over the clean reads
RANK = {"PASS": 0, "WARN": 1, "FAIL": 2}
with open(os.path.join(OUT, "defects_escalation.tsv"), "w") as f:
    f.write("assay\tprofile\tdefect\tgeneric_escalated\tprofile_escalated"
            "\truns\n")
    for a in assays:
        own = OWN.get(a, a)
        for d, m in TARGET.items():
            def esc(p):
                return sum(RANK[by[(a, r, d, p)][m]] >
                           RANK[by[(a, r, "clean", p)][m]] for r in runs[a])
            f.write(f"{a}\t{own}\t{d}\t{esc('generic')}\t{esc(own)}\t"
                    f"{len(runs[a])}\n")

# 3. every profile on every assay's clean runs
profiles = list(dict.fromkeys(p for (_, _, c, p) in by if c == "clean"))
with open(os.path.join(OUT, "cross_profiles.tsv"), "w") as f:
    f.write("assay\tprofile\tmean_flagged\tmean_failed\n")
    for a in assays:
        for p in profiles:
            fl = [sum(map(flagged, by[(a, r, "clean", p)].values()))
                  for r in runs[a]]
            fa = [list(by[(a, r, "clean", p)].values()).count("FAIL")
                  for r in runs[a]]
            f.write(f"{a}\t{p}\t{st.mean(fl):.2f}\t{st.mean(fa):.2f}\n")


# 4. bias extent against the masks (per position, as the grading measures)
def per_base(run):
    rows = list(csv.DictReader(open(os.path.join(OUT, "per_base",
                                                 run + ".tsv")),
                               delimiter="\t"))
    return [(int(r["position"].split("-")[0]), float(r["G"]), float(r["A"]),
             float(r["T"]), float(r["C"])) for r in rows]


with open(os.path.join(OUT, "masks.tsv"), "w") as f:
    f.write("assay\trun\tread_length\tmask_5p\tbias_end_5p_warn"
            "\tbias_end_5p_fail\tmask_3p\tbias_start_3p_warn"
            "\tflagged_beyond_masks\n")
    for a, (m5, m3, meth) in MASK.items():
        for r in runs.get(a, []):
            pb = per_base(r)
            L = len(pb)
            if meth:  # deviation from the mean of the central positions
                core = pb[min(15, L // 4):max(L - 15, 3 * L // 4)]
                mean = [st.mean(x[k] for x in core) for k in range(1, 5)]
                dev = [max(abs(x[k] - mean[k - 1]) for k in range(1, 5))
                       for x in pb]
            else:
                dev = [max(abs(x[2] - x[3]), abs(x[1] - x[4])) for x in pb]
            # last position from the 5' end (within the first 40) above the
            # limit, and the first position from the 3' end (last 40)
            e5w = max([i + 1 for i in range(min(40, L)) if dev[i] > WARN],
                      default=0)
            e5f = max([i + 1 for i in range(min(40, L)) if dev[i] > FAIL],
                      default=0)
            s3w = max([L - i for i in range(max(0, L - 40), L)
                       if dev[i] > WARN], default=0)
            beyond = sum(1 for i in range(m5, L - m3) if dev[i] > WARN)
            f.write(f"{a}\t{r}\t{L}\t{m5}\t{e5w}\t{e5f}\t{m3}\t{s3w}\t"
                    f"{beyond}\n")

# 5. bisulfite/EM conversion annotations
notes = list(csv.DictReader(open(os.path.join(OUT, "annotations.tsv")),
                            delimiter="\t"))
conv = collections.Counter()
for n in notes:
    if n["condition"] != "clean":
        continue
    t = n["annotation"]
    if "conversion" in t.lower() or "depletion" in t.lower():
        key = re.sub(r"[0-9.]+%?", "#", t.split(":", 1)[1].strip())[:90]
        conv[(n["assay"], n["profile"], key)] += 1
with open(os.path.join(OUT, "conversion.tsv"), "w") as f:
    f.write("assay\tprofile\tannotation\truns\n")
    for (a, p, k), c in sorted(conv.items()):
        f.write(f"{a}\t{p}\t{k}\t{c}\n")
print("written to", OUT)
