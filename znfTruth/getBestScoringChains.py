#!/usr/bin/env python

import sys

f = open(sys.argv[1])
s = f.read()
a = map(lambda x: x.split(), s.split("\n"))
bestScoring = {}
for line in a:
    if len(line) < 10:
        continue
    query = line[1]
    score = int(line[9])
    if query in bestScoring:
        if int(bestScoring[query][9]) < score:
            bestScoring[query] = line
    else:
        bestScoring[query] = line

for line in bestScoring.values():
    print " ".join(line)
