#!/usr/bin/env python3
"""Reads whose names change format after the first read: 7-field Illumina
names, then 5-field names with a comment (the tile field of the 7-field
layout falls into the comment, an arbitrary tile number above 32 bits) and
names without fields; negative y coordinates give
negative digits in that number. usage: make_mixed_names.py OUT.fq"""
import random
import sys

rnd = random.Random(3)
with open(sys.argv[1], 'w') as f:
    for k in range(3000):
        v = rnd.random() if k else 0.0
        t = rnd.choice((1101, 1102, 2101))
        x, y = rnd.randrange(1000, 30000), rnd.randrange(-3000, 30000)
        if v < 0.6:
            name = 'A01234:77:HABCDEFXX:1:%d:%d:%d 1:N:0:ACGT' % (t, x, y)
        elif v < 0.9:
            name = 'HWI-ST1:8:%d:%d:%d 1:N:0:ACGT' % (t, x, y)
        else:
            name = 'read%d' % k
        s = ''.join(rnd.choice('ACGT') for _ in range(60))
        q = ''.join(chr(33 + rnd.randrange(2, 41)) for _ in range(60))
        f.write('@%s\n%s\n+\n%s\n' % (name, s, q))
