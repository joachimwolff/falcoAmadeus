#!/usr/bin/env python3
"""Test oracle for the error- and UMI-aware paired duplication.

  aware_oracle.py simulate <seed> <R1.fastq> <R2.fastq>
  aware_oracle.py check <R1.fastq> <R2.fastq|-> <groups.u32> <mismatches>
                  <umi_mismatches> <umi: none|name|r1:S:L|r1r2:S:L>

`simulate` writes read pairs drawn from a few hundred fragments, with
substitutions, N, swapped mates, several read lengths and UMIs (in the first
bases of both reads and in the name). `check` groups the pairs by comparing
every pair with every other one (the definition, without the k-mer blocks)
and compares the groups with those falcoAmadeus wrote; R2 "-" checks
single-end reads.
"""
import random
import sys

import numpy as np


def simulate(seed, out1, out2):
    rnd = random.Random(seed)
    bases = "ACGT"

    def mutate(s, k):
        s = list(s)
        for _ in range(k):
            p = rnd.randrange(len(s))
            s[p] = rnd.choice([b for b in "ACGTN" if b != s[p]])
        return "".join(s)

    pairs = []
    for f in range(300):
        rl = rnd.choice([40, 60, 75, 100, 150])
        a = "".join(rnd.choice(bases) for _ in range(rl))
        b = "".join(rnd.choice(bases) for _ in range(rl))
        u1 = "".join(rnd.choice(bases) for _ in range(6))
        u2 = "".join(rnd.choice(bases) for _ in range(6))
        hu = "".join(rnd.choice(bases) for _ in range(8))
        # low-complexity fragments make large candidate groups
        if f % 50 == 0:
            a = "A" * rl
            b = "AC" * (rl // 2)
        for c in range(rnd.choice([1, 1, 2, 3, 6])):
            x, y = mutate(a, rnd.choice([0, 0, 1, 2, 3, 5])), \
                mutate(b, rnd.choice([0, 0, 1, 2]))
            v1, v2 = mutate(u1, rnd.choice([0, 0, 0, 1])), u2
            w = mutate(hu, rnd.choice([0, 0, 0, 1, 2]))
            x, y = v1 + x[6:], v2 + y[6:]
            if rnd.random() < 0.3:
                x, y = y, x
            pairs.append((x, y, w))
    rnd.shuffle(pairs)
    with open(out1, "w") as o1, open(out2, "w") as o2:
        for i, (x, y, w) in enumerate(pairs):
            name = f"r{i}" if i % 97 == 5 else f"M1:7:FC:1:1:{i}:2:{w}"
            o1.write(f"@{name} 1:N:0\n{x}\n+\n{'I' * len(x)}\n")
            o2.write(f"@{name} 2:N:0\n{y}\n+\n{'I' * len(y)}\n")


def read_fastq(path):
    with open(path) as f:
        lines = f.read().split("\n")
    return [(lines[i][1:], lines[i + 1]) for i in range(0, len(lines) - 3, 4)]


def key(seq, s, l):
    if l:
        seq = seq[:s] + seq[s + l:]
    return seq[:50] if len(seq) > 75 else seq


def dist(x, y):
    if len(x) != len(y):
        return 1 << 20
    return sum(1 for a, b in zip(x, y)
               if a != b or (a == "N") != (b == "N"))


def umi_name(name):
    tok = name.split()[0]
    i = max(tok.rfind(":"), tok.rfind("_"))
    if i < 0:
        return None
    u = tok[i + 1:].replace("+", "").replace("-", "")
    if not u or any(c not in "ACGTNacgtn" for c in u):
        return None
    return u


def check(r1, r2, groups, D, U, umi):
    D, U = int(D), int(U)
    single = r2 == "-"
    m1 = read_fastq(r1)
    m2 = [("", "")] * len(m1) if single else read_fastq(r2)
    s = l = 0
    mode = umi.split(":")[0]
    if mode in ("r1", "r1r2"):
        s, l = int(umi.split(":")[1]) - 1, int(umi.split(":")[2])
    rec = []
    for (n1, a), (_, b) in zip(m1, m2):
        ua = a[s:s + l] if mode in ("r1", "r1r2") else ""
        ub = b[s:s + l] if mode == "r1r2" else ""
        ka = key(a, s, l if mode in ("r1", "r1r2") else 0)
        kb = key(b, s, l if mode == "r1r2" else 0)
        hn = umi_name(n1) if mode == "name" else ""
        rec.append((ka, kb, ua, ub, hn))
    n = len(rec)
    parent = list(range(n))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def udist(x, y):
        if x is None or y is None:
            return 0 if x is y else 1 << 20
        return dist(x, y)

    # mates with keys shorter than 4 bases per block are grouped only when
    # identical (as in falcoAmadeus)
    short = [len(r[0]) < 4 * (D + 1) or
             (not single and len(r[1]) < 4 * (D + 1)) for r in rec]
    for i in range(n):
        ka, kb, ua, ub, hn = rec[i]
        for j in range(i + 1, n):
            kc, kd, uc, ud, hm = rec[j]
            hu = udist(hn, hm)
            if hu > U:
                continue
            Dij, Uij = D, U
            if short[i] or short[j]:
                Dij, Uij = 0, 0
                if hu > 0:
                    continue
            ok = False
            orders = ((kc, kd, uc, ud),) if single else \
                ((kc, kd, uc, ud), (kd, kc, ud, uc))
            for (p, q, up, uq) in orders:
                if dist(ka, p) + dist(kb, q) > Dij:
                    continue
                if hu + udist(ua, up) + udist(ub, uq) > Uij:
                    continue
                ok = True
                break
            if ok:
                ri, rj = find(i), find(j)
                parent[max(ri, rj)] = min(ri, rj)
    want = np.array([find(i) for i in range(n)], dtype=np.uint32)
    # group number = smallest pair index of the group
    first = {}
    for i, g in enumerate(want):
        first.setdefault(int(g), i)
    want = np.array([first[int(g)] for g in want], dtype=np.uint32)
    got = np.fromfile(groups, dtype=np.uint32)
    if len(got) != n or not (got == want).all():
        bad = np.nonzero(got != want)[0][:5] if len(got) == n else []
        print(f"FAIL D={D} U={U} umi={umi}: {len(set(want))} groups expected, "
              f"{len(set(got.tolist()))} found; first differing pairs {bad}")
        return 1
    print(f"ok D={D} U={U} umi={umi}: {len(set(want))} groups of {n} pairs")
    return 0


if __name__ == "__main__":
    if sys.argv[1] == "simulate":
        simulate(int(sys.argv[2]), sys.argv[3], sys.argv[4])
    else:
        sys.exit(check(*sys.argv[2:8]))
