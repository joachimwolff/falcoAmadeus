#!/usr/bin/env python3
"""falcoAmadeus's optical duplicates on the pairs Picard examines.

Picard MarkDuplicates counts optical duplicates among aligned pairs
(both mates mapped, primary); falcoAmadeus among all pairs that would
align. This restricts falcoAmadeus's groups (FA_AWARE_CLUSTERS of the run
with M mismatches) to Picard's pairs, splits them by Picard's orientation
(read 1 reverse or not), and counts optical duplicates by Picard's rule
(same tile, x and y within the distance, x and y in 16 bits, single
linkage). It also reports, for close copies within a falcoAmadeus group,
how many Picard placed at different coordinates and their mapping quality.

usage: compare_on_picard_pairs.py DIR M [DIR is a run_optical_public.sh dir]
"""
import collections
import struct
import sys


def s16(v):
    v &= 0xffff
    return v - 65536 if v >= 32768 else v


def main():
    d, M = sys.argv[1], int(sys.argv[2])
    names = []
    with open(d + '/r1.fq') as f:
        for i, line in enumerate(f):
            if i % 4 == 0:
                names.append(line[1:].split()[0])
    idx = {n: i for i, n in enumerate(names)}
    pos, mq = {}, {}
    for line in open(d + '/picard_pairs.tsv'):
        n, p, q = line.split()
        pos[idx[n]] = p
        mq[idx[n]] = int(q)
    loc = {}
    for i in pos:
        f = names[i].split(':')
        loc[i] = (int(f[4]), s16(int(f[5])), s16(int(f[6])))
    raw = open(d + '/groups_100_%d.bin' % M, 'rb').read()
    grp = struct.unpack('<%dI' % len(names), raw)
    sets = collections.defaultdict(list)
    for i in pos:
        sets[(grp[i], pos[i].split(':')[-1])].append(i)
    gs = collections.Counter(grp[i] for i in pos)
    print('pairs_picard_examines\t%d' % len(pos))
    print('falco_duplicates_on_these\t%d' % (sum(gs.values()) - len(gs)))
    for dist in (100, 2500):
        total, same, diff, diff_mq20 = 0, 0, 0, 0
        for m in sets.values():
            if len(m) < 2:
                continue
            L = sorted(m, key=lambda i: loc[i])
            par = list(range(len(L)))

            def fd(a):
                while par[a] != a:
                    par[a] = par[par[a]]
                    a = par[a]
                return a
            for a in range(len(L)):
                for b in range(a + 1, len(L)):
                    la, lb = loc[L[a]], loc[L[b]]
                    if lb[0] != la[0] or lb[1] - la[1] > dist:
                        break
                    if abs(lb[2] - la[2]) <= dist:
                        ra, rb = fd(a), fd(b)
                        if ra != rb:
                            par[ra] = rb
                        if dist == 100:
                            if pos[L[a]] == pos[L[b]]:
                                same += 1
                            else:
                                diff += 1
                                diff_mq20 += min(mq[L[a]], mq[L[b]]) >= 20
            total += len(L) - len(set(fd(a) for a in range(len(L))))
        print('falco_optical_on_these_%d\t%d' % (dist, total))
        if dist == 100:
            print('close_copies_same_picard_coordinates\t%d' % same)
            print('close_copies_other_picard_coordinates\t%d' % diff)
            print('  of_which_both_mapq_20_or_more\t%d' % diff_mq20)


if __name__ == '__main__':
    main()
