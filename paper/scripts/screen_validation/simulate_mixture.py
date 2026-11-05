#!/usr/bin/env python3
"""Reads drawn from reference sequences in given proportions, for checking
the contamination screen against a known truth.

Usage: simulate_mixture.py <reads> <seed> <out.fastq> name=fraction=fasta ...

Each read (100 bases, either strand, 0.3 % substitutions) comes from a
position drawn uniformly over the named FASTA (concatenated sequences,
windows with N skipped). Writes truth.tsv next to the output: name,
requested fraction, reads drawn.
"""
import gzip
import os
import sys

import numpy as np

N, SEED, OUT = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
L, ERR = 100, 0.003
rng = np.random.default_rng(SEED)
comp = bytes.maketrans(b"ACGTN", b"TGCAN")


def load(path):
    op = gzip.open if path.endswith(".gz") else open
    parts = []
    with op(path, "rb") as f:
        for line in f:
            if not line.startswith(b">"):
                parts.append(line.strip().upper())
    return b"".join(parts)


specs = []
for a in sys.argv[4:]:
    name, frac, path = a.split("=", 2)
    specs.append((name, float(frac), load(path)))
fr = np.array([s[1] for s in specs])
counts = rng.multinomial(N, fr / fr.sum())
letters = np.frombuffer(b"ACGT", dtype=np.uint8)
with open(OUT, "w") as o:
    k = 0
    for (name, _, seq), c in zip(specs, counts):
        n = 0
        while n < c:
            p = int(rng.integers(0, len(seq) - L))
            r = seq[p:p + L]
            if b"N" in r:
                continue
            if rng.random() < 0.5:
                r = r.translate(comp)[::-1]
            a = np.frombuffer(r, dtype=np.uint8).copy()
            m = rng.random(L) < ERR
            a[m] = letters[rng.integers(0, 4, m.sum())]
            o.write(f"@s{k}_{name}\n{a.tobytes().decode()}\n+\n{'I' * L}\n")
            k += 1
            n += 1
with open(os.path.join(os.path.dirname(OUT) or ".", "truth.tsv"), "w") as t:
    t.write("name\tfraction\treads\n")
    for (name, f, _), c in zip(specs, counts):
        t.write(f"{name}\t{f}\t{c}\n")
