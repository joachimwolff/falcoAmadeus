#!/usr/bin/env python3
"""Read pairs without sequencing errors from a random reference, as FASTQ
(R1, R2) and as the matching alignments (SAM), so that copies of a fragment
are both sequence-identical and at the same coordinates: duplicate sets of
falcoAmadeus (sequence) and of Picard MarkDuplicates (alignment) coincide.

Copies of a fragment get Illumina names (instrument:run:flowcell:lane:
tile:x:y, or 5 fields, or an unparseable name) placed near an earlier copy
(within or just beyond a pixel distance, in chains), on another tile, or
anywhere; some copies have the mates swapped (Picard's RF orientation).

usage: make_optical_data.py SEED OUTDIR [MOLECULES]
"""
import random
import sys

COMP = str.maketrans('ACGT', 'TGCA')


def rc(s):
    return s.translate(COMP)[::-1]


def main():
    seed, out = int(sys.argv[1]), sys.argv[2]
    mols = int(sys.argv[3]) if len(sys.argv) > 3 else 4000
    rnd = random.Random(seed)
    G, R = 3_000_000, 100
    ref = ''.join(rnd.choice('ACGT') for _ in range(G))
    with open(out + '/ref.fa', 'w') as f:
        f.write('>chr1\n')
        for i in range(0, G, 80):
            f.write(ref[i:i + 80] + '\n')
    sam = open(out + '/aln.sam', 'w')
    sam.write('@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chr1\tLN:%d\n' % G)
    sam.write('@RG\tID:rg1\tSM:s1\tLB:lib1\n')
    f1, f2 = open(out + '/r1.fq', 'w'), open(out + '/r2.fq', 'w')
    q = 'I' * R
    k = 0
    used, names = set(), set()
    for m in range(mols):
        while True:
            p = rnd.randrange(0, G - 700)
            if p not in used:
                used.add(p)
                break
        L = rnd.randrange(150, 600)
        a, b = ref[p:p + R], ref[p + L - R:p + L]
        r = rnd.random()
        n = 1 if r < 0.55 else (rnd.randrange(2, 8) if r < 0.97
                                else rnd.randrange(50, 2500))
        placed = []
        for c in range(n):
            u = rnd.random()
            if placed and u < 0.45:  # near an earlier copy
                t, x, y = rnd.choice(placed)
                dx = rnd.choice((rnd.randrange(0, 101), rnd.randrange(95, 140),
                                 rnd.randrange(2000, 2600)))
                dy = rnd.choice((rnd.randrange(0, 101), rnd.randrange(95, 140),
                                 rnd.randrange(0, 2600)))
                x, y = x + rnd.choice((-1, 1)) * dx, y + rnd.choice((-1, 1)) * dy
            elif placed and u < 0.6:  # same spot, another tile
                t, x, y = rnd.choice(placed)
                t = t + rnd.choice((1, 100, 1000))
            else:
                t = rnd.choice((1101, 1102, 1201, 2101, 2214))
                x, y = rnd.randrange(1000, 30000), rnd.randrange(1000, 30000)
            while (t, x, y) in names:  # unique names
                x += 1
            names.add((t, x, y))
            placed.append((t, x, y))
            v = rnd.random()
            if v < 0.8:
                name = 'A01234:77:HABCDEFXX:1:%d:%d:%d' % (t, x, y)
            elif v < 0.95:
                name = 'HWI-ST1:8:%d:%d:%d' % (t, x, y)
            else:
                name = 'read%d' % k
            name += '_%d' % k if v >= 0.95 else ''
            k += 1
            swap = rnd.random() < 0.3
            s1, s2 = (a, rc(b)) if not swap else (rc(b), a)
            # FR: read 1 forward at p, read 2 reverse at p+L-R
            # RF: read 1 reverse at p+L-R, read 2 forward at p
            f1.write('@%s 1:N:0:ACGT\n%s\n+\n%s\n' % (name, s1, q))
            f2.write('@%s 2:N:0:ACGT\n%s\n+\n%s\n' % (name, s2, q))
            pf, pr = p + 1, p + L - R + 1
            tl = L
            if not swap:
                sam.write('%s\t99\tchr1\t%d\t60\t%dM\t=\t%d\t%d\t%s\t%s\tRG:Z:rg1\n'
                          % (name, pf, R, pr, tl, a, q))
                sam.write('%s\t147\tchr1\t%d\t60\t%dM\t=\t%d\t%d\t%s\t%s\tRG:Z:rg1\n'
                          % (name, pr, R, pf, -tl, b, q))
            else:
                sam.write('%s\t83\tchr1\t%d\t60\t%dM\t=\t%d\t%d\t%s\t%s\tRG:Z:rg1\n'
                          % (name, pr, R, pf, -tl, b, q))
                sam.write('%s\t163\tchr1\t%d\t60\t%dM\t=\t%d\t%d\t%s\t%s\tRG:Z:rg1\n'
                          % (name, pf, R, pr, tl, a, q))
    print(k)


if __name__ == '__main__':
    main()
