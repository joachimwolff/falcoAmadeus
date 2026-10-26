#!/usr/bin/env python3
"""Sequence duplication of read pairs against Picard MarkDuplicates on the
same pairs.

Usage: compare_picard.py <marked.bam> <R1.fastq> <R2.fastq> <out_dir>

Compared pairs: both mates mapped as primary alignments, properly paired,
mapping quality >= 30 for both. A pair's sequence key is the FastQC
duplication key of each mate, in either mate order (the first 50 bases of reads
longer than 75 bases, otherwise the whole read), as in falcoAmadeus. Writes
picard_comparison.tsv and the compared pairs as compared_R1/R2.fastq.
"""
import hashlib
import os
import sys

import numpy as np
import pysam

BAM, R1, R2, OUT = sys.argv[1:5]
MIN_MAPQ = 30


def read_fastq(path):
    with open(path, "rb") as f:
        while True:
            h = f.readline()
            if not h:
                return
            s = f.readline().rstrip(b"\n")
            p = f.readline()
            q = f.readline()
            yield h, s, p, q


def key(seq):
    return seq[:50] if len(seq) > 75 else seq


# pass 1: pair keys in file order
keys = []
for (h1, s1, _, _), (h2, s2, _, _) in zip(read_fastq(R1), read_fastq(R2)):
    a, b = sorted((key(s1), key(s2)))  # either mate order
    d = hashlib.blake2b(a + b"|" + b, digest_size=8).digest()
    keys.append(int.from_bytes(d, "little"))
keys = np.array(keys, dtype=np.uint64)
n = len(keys)

# pass 2: alignment state per mate (read index from the name SRR.N)
mapq = np.zeros((2, n), dtype=np.uint8)
mapped = np.zeros((2, n), dtype=bool)
proper = np.zeros((2, n), dtype=bool)
dup = np.zeros((2, n), dtype=bool)
with pysam.AlignmentFile(BAM, "rb") as bam:
    for r in bam.fetch(until_eof=True):
        if r.is_secondary or r.is_supplementary:
            continue
        i = int(r.query_name.rsplit(".", 1)[1]) - 1
        m = 0 if r.is_read1 else 1
        mapped[m, i] = not r.is_unmapped
        mapq[m, i] = r.mapping_quality
        proper[m, i] = r.is_proper_pair
        dup[m, i] = r.is_duplicate

S = (mapped.all(0) & proper.all(0) & (mapq >= MIN_MAPQ).all(0))
ns = int(S.sum())
picard_pair = dup.all(0)
picard_dup_s = int((picard_pair & S).sum())
mate_disagree = int(((dup[0] != dup[1]) & S).sum())

ks = keys[S]
uniq, inv, counts = np.unique(ks, return_inverse=True, return_counts=True)
seq_extra = ns - len(uniq)  # pairs that repeat an earlier pair's sequences
in_group = counts[inv] > 1
pf = picard_pair[S]
# Picard-flagged pairs that have a sequence-identical partner
flagged_with_partner = int((pf & in_group).sum())
# extra copies of sequence groups that Picard also flagged (per group: at
# most size - 1)
flagged_per_group = np.bincount(inv, weights=pf.astype(np.int64),
                                minlength=len(uniq))
extra_flagged = int(np.minimum(flagged_per_group, counts - 1).sum())

rows = [
    ("pairs_total", n),
    ("pairs_compared", ns),
    ("min_mapq", MIN_MAPQ),
    ("picard_duplicate_pairs", picard_dup_s),
    ("picard_duplicate_rate_pct", 100.0 * picard_dup_s / ns),
    ("picard_mates_disagree", mate_disagree),
    ("sequence_duplicate_pairs", seq_extra),
    ("sequence_duplicate_rate_pct", 100.0 * seq_extra / ns),
    ("picard_flagged_with_sequence_partner", flagged_with_partner),
    ("picard_flagged_with_sequence_partner_pct",
     100.0 * flagged_with_partner / max(1, picard_dup_s)),
    ("sequence_extra_copies_flagged_by_picard", extra_flagged),
    ("sequence_extra_copies_flagged_by_picard_pct",
     100.0 * extra_flagged / max(1, seq_extra)),
]
with open(os.path.join(OUT, "picard_comparison.tsv"), "w") as f:
    f.write("metric\tvalue\n")
    for k, v in rows:
        f.write(f"{k}\t{v:.4f}\n" if isinstance(v, float) else f"{k}\t{v}\n")

# the compared pairs as FASTQ, in file order (for falcoAmadeus)
with open(os.path.join(OUT, "compared_R1.fastq"), "wb") as o1, \
     open(os.path.join(OUT, "compared_R2.fastq"), "wb") as o2:
    for i, (a, b) in enumerate(zip(read_fastq(R1), read_fastq(R2))):
        if S[i]:
            o1.write(a[0] + a[1] + b"\n" + a[2] + a[3])
            o2.write(b[0] + b[1] + b"\n" + b[2] + b[3])
for k, v in rows:
    print(f"{k}\t{v}")
