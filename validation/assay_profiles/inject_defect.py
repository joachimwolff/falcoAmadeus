#!/usr/bin/env python3
"""Inject one controlled QC defect into a FASTQ file (deterministic).

Defects, each applied to a fraction of reads chosen with a seeded RNG:

  adapter     overwrite the 3' end with the Illumina TruSeq adapter
              (AGATCGGAAGAGC...) starting at 60% of the read length
  qualtail    set the last 20 quality values to Q2 ('#')
  midbias     overwrite read positions 20-29 with 'G'; a composition defect
              downstream of every 5' mask used by the assay profiles
  gccontam    replace the read with a random sequence of 75% GC
              (simulates contamination by a foreign high-GC source)
  dupspike    replace the read with one of 20 template reads, so that the
              chosen fraction of the file is high-copy duplicates
  ncontent    set positions 10-19 to 'N'

Usage: inject_defect.py IN.fastq OUT.fastq DEFECT FRACTION [--seed S]
"""
import argparse
import random

ADAPTER = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTGCTTG"


def gc_seq(rng, n):
    return "".join(rng.choice("GC") if rng.random() < 0.75 else rng.choice("AT")
                   for _ in range(n))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("inp")
    ap.add_argument("out")
    ap.add_argument("defect", choices=["adapter", "qualtail", "midbias",
                                       "gccontam", "dupspike", "ncontent"])
    ap.add_argument("fraction", type=float)
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    rng = random.Random(a.seed)
    templates = []

    with open(a.inp) as fi, open(a.out, "w") as fo:
        while True:
            h = fi.readline()
            if not h:
                break
            s = fi.readline().rstrip("\n")
            p = fi.readline()
            q = fi.readline().rstrip("\n")
            n = len(s)
            if a.defect == "dupspike" and len(templates) < 20 and n >= 30:
                templates.append((s, q))
            if n >= 30 and rng.random() < a.fraction:
                if a.defect == "adapter":
                    k = int(0.6 * n)
                    s = (s[:k] + ADAPTER * (1 + n // len(ADAPTER)))[:n]
                elif a.defect == "qualtail":
                    q = q[:n - 20] + "#" * 20
                elif a.defect == "midbias":
                    s = s[:20] + "G" * 10 + s[30:]
                elif a.defect == "gccontam":
                    s = gc_seq(rng, n)
                elif a.defect == "dupspike" and templates:
                    s, q = templates[rng.randrange(len(templates))]
                elif a.defect == "ncontent":
                    s = s[:10] + "N" * 10 + s[20:]
            fo.write(f"{h}{s}\n{p}{q}\n")


if __name__ == "__main__":
    main()
