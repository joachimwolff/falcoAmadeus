#!/usr/bin/env python3
"""Groups of copies from falcoAmadeus (error- and UMI-aware paired mode)
against Picard MarkDuplicates on the same pairs.

Usage: compare_picard_groups.py <marked.bam> <groups.u32> <label> <out.tsv>

<groups.u32> holds one 32-bit group number per compared pair, in the order
of compared_R1/R2.fastq (written by falcoAmadeus with FA_AWARE_CLUSTERS set,
run on the compared pairs). The compared pairs are selected as in
compare_picard.py. Appends one row to <out.tsv>.
"""
import os
import sys

import numpy as np
import pysam

BAM, GROUPS, LABEL, OUT = sys.argv[1:5]
MIN_MAPQ = 30

with pysam.AlignmentFile(BAM, "rb") as bam:
    n = 0
    rows = []
    for r in bam.fetch(until_eof=True):
        if r.is_secondary or r.is_supplementary:
            continue
        i = int(r.query_name.rsplit(".", 1)[1]) - 1
        rows.append((i, 0 if r.is_read1 else 1, not r.is_unmapped,
                     r.mapping_quality, r.is_proper_pair, r.is_duplicate))
        n = max(n, i + 1)
mapq = np.zeros((2, n), dtype=np.uint8)
mapped = np.zeros((2, n), dtype=bool)
proper = np.zeros((2, n), dtype=bool)
dup = np.zeros((2, n), dtype=bool)
for i, m, mp, q, pp, d in rows:
    mapped[m, i], mapq[m, i], proper[m, i], dup[m, i] = mp, q, pp, d
S = mapped.all(0) & proper.all(0) & (mapq >= MIN_MAPQ).all(0)
pf = dup.all(0)[S]
g = np.fromfile(GROUPS, dtype=np.uint32)
assert len(g) == len(pf), (len(g), len(pf))

uniq, inv, counts = np.unique(g, return_inverse=True, return_counts=True)
extra = len(g) - len(uniq)
flagged = np.bincount(inv, weights=pf.astype(np.int64), minlength=len(uniq))
extra_flagged = int(np.minimum(flagged, counts - 1).sum())
picard = int(pf.sum())
row = {
    "label": LABEL,
    "pairs_compared": len(g),
    "picard_duplicate_rate_pct": 100.0 * picard / len(g),
    "falcoamadeus_duplicate_rate_pct": 100.0 * extra / len(g),
    "falcoamadeus_extra_copies": extra,
    "extra_copies_flagged_by_picard_pct": 100.0 * extra_flagged / max(1, extra),
    "picard_duplicates_found_pct": 100.0 * extra_flagged / max(1, picard),
}
new = not os.path.exists(OUT)
with open(OUT, "a") as f:
    if new:
        f.write("\t".join(row) + "\n")
    f.write("\t".join(f"{v:.4f}" if isinstance(v, float) else str(v)
                      for v in row.values()) + "\n")
print(row)
