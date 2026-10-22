#!/usr/bin/env python3
"""Deterministic paired-end FASTQ simulator with a known duplicate fraction.

Writes R1/R2 FASTQ files of N read pairs in which exactly round(N * dup) pairs
are copies of an earlier distinct pair (both mates identical), so the true
sequence-pair duplicate rate is dup. The placement of the copies is
controlled to probe order dependence:

  late      all distinct pairs first, all copies at the end
  early     copies follow their source pair at the start of the file
  shuffled  uniformly random order (seeded)

Additional modes:
  --shared-r1    copies repeat R1 only; R2 is always new (true paired rate 0)
  --reorder-r2   shuffle the R2 records (mates no longer in matching order)

Output is fully determined by the arguments and --seed.
"""
import argparse
import random


def rand_seq(rng, n):
    return "".join(rng.choice("ACGT") for _ in range(n))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("--pairs", type=int, required=True)
    ap.add_argument("--dup", type=float, required=True,
                    help="fraction of pairs that are copies (0..1)")
    ap.add_argument("--placement", choices=["late", "early", "shuffled"],
                    default="shuffled")
    ap.add_argument("--length", type=int, default=100)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--shared-r1", action="store_true")
    ap.add_argument("--reorder-r2", action="store_true")
    ap.add_argument("--prefix", required=True,
                    help="writes <prefix>_R1.fastq and <prefix>_R2.fastq")
    a = ap.parse_args()

    rng = random.Random(a.seed)
    n_copies = round(a.pairs * a.dup)
    n_distinct = a.pairs - n_copies
    distinct = [(rand_seq(rng, a.length), rand_seq(rng, a.length))
                for _ in range(n_distinct)]
    sources = [rng.randrange(n_distinct) for _ in range(n_copies)]

    def copy_of(src):
        r1, r2 = distinct[src]
        return (r1, rand_seq(rng, a.length)) if a.shared_r1 else (r1, r2)

    if a.placement == "late":
        pairs = distinct + [copy_of(s) for s in sources]
    elif a.placement == "early":
        by_src = {}
        for s in sources:
            by_src.setdefault(s, []).append(copy_of(s))
        # sources with copies first, each immediately followed by its copies
        head, tail = [], []
        for i, p in enumerate(distinct):
            if i in by_src:
                head.append(p)
                head.extend(by_src[i])
            else:
                tail.append(p)
        pairs = head + tail
    else:
        pairs = distinct + [copy_of(s) for s in sources]
        rng.shuffle(pairs)

    r2_order = list(range(len(pairs)))
    if a.reorder_r2:
        rng.shuffle(r2_order)

    qual = "I" * a.length
    with open(a.prefix + "_R1.fastq", "w") as f1, \
         open(a.prefix + "_R2.fastq", "w") as f2:
        for i, (r1, _) in enumerate(pairs):
            f1.write(f"@sim.{i} 1\n{r1}\n+\n{qual}\n")
        for i in r2_order:
            f2.write(f"@sim.{i} 2\n{pairs[i][1]}\n+\n{qual}\n")


if __name__ == "__main__":
    main()
