#!/usr/bin/env python3
"""Trimming preview against cutadapt.

For each FASTQ file: cutadapt -q C (C = 10, 20, 30) is run on the reads the
preview examined; the bases removed, reads trimmed and reads shorter than 20,
36 and 50 bases are counted from cutadapt's output and formatted as the
preview prints them (%g). The quality rows must be equal. The adapter row is
an estimate; it is compared with cutadapt -a ADAPTER -q 20 and the
difference is reported.

usage: compare_cutadapt.py FALCOAMADEUS CUTADAPT OUTDIR FASTQ...
"""
import os
import subprocess
import sys


def records(path):
    with open(path) as f:
        while True:
            h = f.readline()
            if not h:
                return
            s = f.readline().rstrip('\n')
            f.readline()
            q = f.readline().rstrip('\n')
            yield s, q


def counts(lengths_before, path):
    removed = trimmed = 0
    short = [0, 0, 0]
    bases = 0
    for (n, (s, _)) in zip(lengths_before, records(path)):
        k = len(s)
        bases += n
        removed += n - k
        trimmed += k < n
        for i, m in enumerate((20, 36, 50)):
            short[i] += k < m
    return bases, removed, trimmed, short


def pct(x, n):
    return '%g' % (100.0 * x / n if n else 0.0)


def preview_rows(path):
    rows, adapter = {}, {}
    with open(path) as f:
        for line in f:
            f_ = line.rstrip('\n').split('\t')
            if f_[0] in ('10', '20', '30') and len(f_) == 6:
                rows[f_[0]] = f_[1:]
            elif f_[0] == 'Adapter':
                adapter['name'] = f_[1:]
            elif len(f_) == 2 and 'with quality 20' in f_[0]:
                adapter[f_[0]] = f_[1]
    return rows, adapter


def main():
    fa, cutadapt, out = sys.argv[1], sys.argv[2], sys.argv[3]
    os.makedirs(out, exist_ok=True)
    ok = True
    print('file\tcutoff\tfield\tpreview\tcutadapt\tequal')
    for fq in sys.argv[4:]:
        name = os.path.basename(fq)
        d = os.path.join(out, name + '.fa')
        os.makedirs(d, exist_ok=True)
        subprocess.run([fa, '-q', '--trim-preview', '--trim-preview-reads',
                        '0', '-o', d, fq], check=True)
        rows, adapter = preview_rows(os.path.join(d, 'trimming_preview.txt'))
        lens = [len(s) for s, _ in records(fq)]
        n = len(lens)
        for c in ('10', '20', '30'):
            t = os.path.join(out, name + '.q' + c + '.fastq')
            subprocess.run([cutadapt, '-q', c, '-o', t, fq], check=True,
                           stdout=subprocess.DEVNULL,
                           stderr=subprocess.DEVNULL)
            bases, removed, trimmed, short = counts(lens, t)
            ref = [pct(removed, bases), pct(trimmed, n)] + \
                [pct(x, n) for x in short]
            for field, a, b in zip(('bases_removed', 'reads_trimmed',
                                    'short20', 'short36', 'short50'),
                                   rows[c], ref):
                eq = a == b
                ok &= eq
                print('%s\t%s\t%s\t%s\t%s\t%s' % (name, c, field, a, b, eq))
        if len(adapter.get('name', [])) == 2:
            seq = adapter['name'][1]
            t = os.path.join(out, name + '.a.fastq')
            subprocess.run([cutadapt, '-a', seq, '-q', '20', '-o', t, fq],
                           check=True, stdout=subprocess.DEVNULL,
                           stderr=subprocess.DEVNULL)
            bases, removed, trimmed, short = counts(lens, t)
            ref = {'Bases removed with quality 20 (%)': pct(removed, bases),
                   'Reads trimmed with quality 20 (%)': pct(trimmed, n),
                   'Shorter than 36 with quality 20 (%)': pct(short[1], n)}
            for k, v in ref.items():
                print('%s\tadapter+20\t%s\t%s\t%s\testimate' %
                      (name, k, adapter.get(k), v))
    sys.exit(0 if ok else 1)


if __name__ == '__main__':
    main()
