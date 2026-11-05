#!/usr/bin/env python3
"""Independent computation of umi_qc.txt (without the File line).

usage: umi_oracle.py CLUSTERS R1 [R2] --umi-r1 S:L --umi-r2 S:L --umi-in-name

CLUSTERS: the group of each read (pair) as 32-bit integers
(FA_AWARE_CLUSTERS of the same run); it gives the families.
"""
import argparse
import collections
import struct

CODE = {'A': 0, 'C': 1, 'G': 2, 'T': 3}


def g(x):
    return '%g' % x


def pct(a, b):
    return g(100.0 * a / b if b else 0.0)


def records(path):
    with open(path) as f:
        while True:
            h = f.readline()
            if not h:
                return
            s = f.readline().rstrip('\n').rstrip('\r')
            f.readline()
            q = f.readline().rstrip('\n').rstrip('\r')
            yield h.rstrip('\n'), s, q


def name_umi(h):
    b = 1 if h.startswith('@') else 0
    e = b
    while e < len(h) and h[e] not in ' \t\r':
        e += 1
    if e - b >= 2 and h[e - 2] in '/.' and h[e - 1] in '12':
        e -= 2
    f = e
    while f > b and h[f - 1] not in ':_':
        f -= 1
    if f == b:
        return None
    out = ''
    for c in h[f:e]:
        if c in '+-':
            continue
        if c.upper() not in 'ACGTN':
            return None
        out += c
    return out if out else None


def norm(u):
    return ''.join(c.upper() if c.upper() in CODE else 'N' for c in u)


def packed(u):
    return sum(CODE[c] << (2 * p) for p, c in enumerate(u))


def seq_section(name, umis, unit):
    """umis: per read the UMI string ('' = none present)"""
    out = ['#' + name, '%s %s\t%d' % (name, unit, len(umis))]
    present = [u for u in umis if u]
    out.append('%s present\t%d\t%s' % (name, len(present),
                                       pct(len(present), len(umis))))
    if not present:
        return out, 0
    lens = [len(u) for u in present]
    lo, hi = min(lens), max(lens)
    out.append('%s length\t%s' % (name, str(lo) if lo == hi else
                                  '%d-%d' % (lo, hi)))
    withn = [u for u in present if 'N' in u]
    clean = [u for u in present if 'N' not in u]
    homo = [u for u in clean if len(u) >= 2 and len(set(u)) == 1]
    out.append('%s with N\t%d\t%s' % (name, len(withn),
                                      pct(len(withn), len(present))))
    out.append('%s homopolymer\t%d\t%s' % (name, len(homo),
                                           pct(len(homo), len(present))))
    cnt = collections.Counter(clean)
    out.append('%s distinct\t%d' % (name, len(cnt)))
    if lo == hi and hi <= 31:
        out.append('%s share of the possible UMIs used (%%)\t%s' %
                   (name, g(100.0 * len(cnt) / 4 ** hi)))
    if cnt:
        best = max(cnt.values())
        if lo == hi:
            top = min((u for u in cnt if cnt[u] == best), key=packed)
        else:
            top = min((u for u in cnt if cnt[u] == best),
                      key=lambda u: (packed(u), len(u)))
    else:
        top, best = '', 0
    out.append('%s most frequent\t%s\t%s' % (name, top,
                                             pct(best, len(clean))))
    out.append('#%s base composition (%%)' % name)
    out.append('#Position\tA\tC\tG\tT\tN')
    for p in range(hi):
        c = collections.Counter(u[p] for u in present if len(u) > p)
        t = sum(c.values())
        out.append('\t'.join([str(p + 1)] + [pct(c[b], t) for b in 'ACGTN']))
    return out, len(cnt)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('clusters')
    ap.add_argument('r1')
    ap.add_argument('r2', nargs='?')
    ap.add_argument('--umi-r1')
    ap.add_argument('--umi-r2')
    ap.add_argument('--umi-in-name', action='store_true')
    a = ap.parse_args()
    reg = lambda r: (int(r.split(':')[0]) - 1, int(r.split(':')[1]))
    unit = 'pairs' if a.r2 else 'reads'
    recs1 = list(records(a.r1))
    recs2 = list(records(a.r2)) if a.r2 else []
    quals1 = [q for _, _, q in recs1]
    allq = ''.join(quals1 + [q for _, _, q in recs2])
    base = 33 if min(map(ord, allq)) < 64 else 64

    def inread(recs, r):
        s, l = reg(r)
        return [norm(seq[s:s + l]) for _, seq, _ in recs]

    out, tags = [], None
    if a.umi_in_name:
        u = [name_umi(h) or '' for h, _, _ in recs1]
        o, d = seq_section('Read name UMI', [norm(x) for x in u], unit)
        out += o
        tags = d
    u1 = inread(recs1, a.umi_r1) if a.umi_r1 else None
    u2 = inread(recs2, a.umi_r2) if (a.umi_r2 and a.r2) else None
    if u1 is not None:
        o, d = seq_section('Read 1 UMI', u1, unit)
        out += o
        tags = d if tags is None else tags
    if u2 is not None:
        o, d = seq_section('Read 2 UMI', u2, unit)
        out += o
        tags = d if tags is None else tags
    if u1 is not None and u2 is not None and not a.umi_in_name:
        tags = len({(x, y) for x, y in zip(u1, u2)
                    if x and y and 'N' not in x and 'N' not in y})
    data = open(a.clusters, 'rb').read()
    grp = struct.unpack('<%dI' % (len(data) // 4), data)
    sizes = collections.Counter(collections.Counter(grp).values())
    fam = sum(sizes.values())
    mem = sum(k * v for k, v in sizes.items())
    ge = lambda m: sum(v for k, v in sizes.items() if k >= m)
    rge = lambda m: sum(k * v for k, v in sizes.items() if k >= m)
    seen, median = 0, 0
    for k in sorted(sizes):
        seen += sizes[k]
        if 2 * seen >= fam:
            median = k
            break
    out += ['#Families (groups of copies: same UMI and sequence within the '
            'allowed mismatches)',
            'Families\t%d' % fam,
            'Mean family size\t%s' % g(mem / fam if fam else 0.0),
            'Median family size\t%d' % median,
            'Largest family\t%d' % (max(sizes) if sizes else 0),
            'Families of 2 or more\t%d\t%s' % (ge(2), pct(ge(2), fam)),
            'Families of 3 or more\t%d\t%s' % (ge(3), pct(ge(3), fam)),
            'Reads in families of 2 or more\t%d\t%s' % (rge(2),
                                                       pct(rge(2), mem)),
            'Reads in families of 3 or more\t%d\t%s' % (rge(3),
                                                       pct(rge(3), mem)),
            'Distinct molecule tags\t%d' % (tags or 0),
            'Families per distinct molecule tag\t%s' %
            g(fam / tags if tags else 0.0),
            '#Family size (from)\tFamilies\tPercentage of families\t'
            'Percentage of %s' % unit]
    frm = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 21, 51, 101, 1001]
    bf, br = [0] * len(frm), [0] * len(frm)
    for k, v in sizes.items():
        b = max(i for i, f in enumerate(frm) if k >= f)
        bf[b] += v
        br[b] += k * v
    for i, f in enumerate(frm):
        out.append('%d\t%d\t%s\t%s' % (f, bf[i], pct(bf[i], fam),
                                       pct(br[i], mem)))
    for name, recs, r in (('Read 1 UMI', recs1, a.umi_r1),
                          ('Read 2 UMI', recs2, a.umi_r2 if a.r2 else None)):
        if not r:
            continue
        s, l = reg(r)
        qs = [q[s:s + l] for _, _, q in recs[:1000000] if len(q) > s]
        low = sum(1 for q in qs if min(ord(c) & 127 for c in q) - base < 20)
        out += ['#%s base quality (the first reads of the file)' % name,
                '%s qualities examined\t%d' % (name, len(qs)),
                '%s with a base below Phred 20\t%d\t%s' % (name, low,
                                                          pct(low, len(qs))),
                '#Position\tMean quality\tBelow Phred 20 (%)']
        for p in range(l):
            v = [(ord(q[p]) & 127) - base for q in qs if len(q) > p]
            out.append('%d\t%s\t%s' % (p + 1, g(sum(v) / len(v) if v else 0.0),
                                       pct(sum(1 for x in v if x < 20),
                                           len(v))))
    print('\n'.join(out))


if __name__ == '__main__':
    main()
