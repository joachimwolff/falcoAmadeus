#!/usr/bin/env python3
"""Sequence differences within Picard duplicate sets.

Usage: picard_set_mismatches.py <marked_with_sets.bam> <R1.fastq> <R2.fastq>
                                <out.tsv>

The BAM comes from MarkDuplicates with TAG_DUPLICATE_SET_MEMBERS=true (DI
tag: duplicate set). For every pair of a set (both mates primary, properly
paired, mapping quality >= 30) the FastQC keys of both mates (first 50
bases) are compared with those of the set's first pair in file order; the
table counts the mismatching bases (read 1 + read 2), also allowing the
mates to be swapped (the other strand of the same fragment), and the pairs
whose keys contain an N.
"""
import collections
import sys

import pysam

BAM, R1, R2, OUT = sys.argv[1:5]
MIN_MAPQ = 30

sets = collections.defaultdict(dict)  # DI -> {read index: ok flags}
state = {}
with pysam.AlignmentFile(BAM, "rb") as bam:
    for r in bam.fetch(until_eof=True):
        if r.is_secondary or r.is_supplementary or not r.has_tag("DI"):
            continue
        i = int(r.query_name.rsplit(".", 1)[1]) - 1
        ok = (not r.is_unmapped and r.is_proper_pair
              and r.mapping_quality >= MIN_MAPQ)
        st = state.setdefault(i, [r.get_tag("DI"), True, 0])
        st[1] = st[1] and ok
        st[2] += 1
members = collections.defaultdict(list)
for i, (di, ok, n) in state.items():
    if ok and n == 2:
        members[di].append(i)
wanted = {i for m in members.values() if len(m) > 1 for i in m}


def fq(path):
    with open(path) as f:
        for k, line in enumerate(f):
            if k % 4 == 1:
                yield line.rstrip("\n")


seq = {}
for i, (a, b) in enumerate(zip(fq(R1), fq(R2))):
    if i in wanted:
        seq[i] = (a[:50] if len(a) > 75 else a, b[:50] if len(b) > 75 else b)

hist = collections.Counter()
swap = collections.Counter()  # either mate order
swapped = 0
with_n = 0
for m in members.values():
    if len(m) < 2:
        continue
    m.sort()
    r1, r2 = seq[m[0]]
    for i in m[1:]:
        a, b = seq[i]
        if "N" in a + b + r1 + r2:
            with_n += 1
        if len(a) != len(r1) or len(b) != len(r2):
            hist["length differs"] += 1
            continue
        d = sum(x != y for x, y in zip(a, r1)) + sum(x != y for x, y in zip(b, r2))
        hist[d if d <= 5 else ">5"] += 1
        # the same fragment read from the other strand: mates swapped
        if len(a) == len(r2) and len(b) == len(r1):
            ds = sum(x != y for x, y in zip(a, r2)) + sum(x != y for x, y in zip(b, r1))
            swap[min(d, ds) if min(d, ds) <= 5 else ">5"] += 1
            if ds < d:
                swapped += 1
total = sum(hist.values())
with open(OUT, "w") as f:
    f.write("mismatches\tpairs\tpct\tpairs_either_order\tpct_either_order\n")
    for k in [0, 1, 2, 3, 4, 5, ">5", "length differs"]:
        f.write(f"{k}\t{hist[k]}\t{100.0 * hist[k] / max(1, total):.2f}"
                f"\t{swap[k]}\t{100.0 * swap[k] / max(1, total):.2f}\n")
    f.write(f"mates_swapped\t{swapped}\t{100.0 * swapped / max(1, total):.2f}\n")
    f.write(f"with_N\t{with_n}\t{100.0 * with_n / max(1, total):.2f}\n")
print(open(OUT).read())
