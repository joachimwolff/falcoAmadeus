#!/usr/bin/env python3
"""Paired-end library with known duplication, for comparing duplicate
counters.

Usage: simulate_duplicates.py <pairs> <seed> <R1.fastq> <R2.fastq> <truth.tsv>

Molecules are fragments (200 to 500 bp) of a random 50 Mb genome. Each end
of a molecule carries its own 8 nt UMI in the adapter, read as the first 8
bases of the read that starts at that end (92 genomic bases follow). A tenth
of the molecules start and end where another molecule does (the same
fragment from two molecules: natural duplicates that only the UMIs
separate). Every molecule is sequenced 1 + Poisson(lam) times (its PCR
copies), lam chosen for the requested number of pairs; each copy is read
from either strand, so read 1 starts at either end (mates and UMIs swap),
and carries substitution errors at 0.5 % per base, in the UMIs as well.
Reads are 100 bp, shuffled. truth.tsv lists, for pairs, the PCR copies
(same molecule) and the copies of a fragment whatever the molecule, and for
read 1 alone the same for reads that start at the same end (a read 1 shows
one end of the fragment only).
"""
import sys

import numpy as np

PAIRS, SEED = int(sys.argv[1]), int(sys.argv[2])
OUT1, OUT2, TRUTH = sys.argv[3:6]
G, L, U = 50_000_000, 100, 8
ERR = 0.005
rng = np.random.default_rng(SEED)
genome = rng.integers(0, 4, G, dtype=np.uint8)
comp = np.array([3, 2, 1, 0], dtype=np.uint8)
letters = np.frombuffer(b"ACGT", dtype=np.uint8)

lam = 0.43
mol = int(PAIRS / (1 + lam))
copies = 1 + rng.poisson(lam, mol)
# trim or pad the copy numbers to exactly PAIRS reads
diff = copies.sum() - PAIRS
while diff > 0:
    i = rng.integers(0, mol)
    if copies[i] > 1:
        copies[i] -= 1
        diff -= 1
if diff < 0:
    copies[rng.integers(0, mol, -diff)] += 1
flen = rng.integers(200, 501, mol)
start = rng.integers(0, G - 501, mol)
shared = rng.random(mol) < 0.10  # reuse an earlier molecule's fragment
src = rng.integers(0, mol, mol)
start[shared] = start[src[shared]]
flen[shared] = flen[src[shared]]
umi = rng.integers(0, 4, (mol, 2, U), dtype=np.uint8)  # one per end

reads = np.repeat(np.arange(mol), copies)
rng.shuffle(reads)


def mutate(a):
    m = rng.random(a.shape) < ERR
    a = a.copy()
    a[m] = (a[m] + rng.integers(1, 4, m.sum(), dtype=np.uint8)) % 4
    return a


q = b"I" * L
G_L = L - U  # genomic bases per read
strand = np.zeros(PAIRS, dtype=np.uint8)
with open(OUT1, "wb") as o1, open(OUT2, "wb") as o2:
    for k, m in enumerate(reads):
        s, f = start[m], flen[m]
        fwd = np.concatenate([umi[m, 0], genome[s:s + G_L]])
        rev = np.concatenate([umi[m, 1], comp[genome[s + f - G_L:s + f][::-1]]])
        a, b = mutate(fwd), mutate(rev)
        if rng.random() < 0.5:  # read 1 starts at the other end
            a, b = b, a
            strand[k] = 1
        n = f"@p.{k + 1}\n".encode()
        o1.write(n + letters[a].tobytes() + b"\n+\n" + q + b"\n")
        o2.write(n + letters[b].tobytes() + b"\n+\n" + q + b"\n")

frag = (start * 1000 + flen).astype(np.int64)


def rate(keys):
    return 100 * (PAIRS - len(np.unique(keys))) / PAIRS


with open(TRUTH, "w") as f:
    f.write("metric\tvalue\n")
    f.write(f"pairs\t{PAIRS}\n")
    f.write(f"molecules\t{mol}\n")
    f.write(f"pe_pcr_duplicate_rate_pct\t{rate(reads):.4f}\n")
    f.write(f"pe_fragment_duplicate_rate_pct\t{rate(frag[reads]):.4f}\n")
    # read 1 alone: the same molecule (or fragment) read from the same end
    f.write(f"se_pcr_duplicate_rate_pct\t{rate(reads.astype(np.int64) * 2 + strand):.4f}\n")
    # a read 1 shows only the position and strand where it starts
    r1start = np.where(strand == 0, start[reads], start[reads] + flen[reads])
    f.write(f"se_start_duplicate_rate_pct\t{rate(r1start.astype(np.int64) * 2 + strand):.4f}\n")
