#!/usr/bin/env python
# Usage: pastaIdsToOriginalNames.py fastaFile renameFile
import sys
from sonLib.bioio import system, fastaRead, fastaWrite

fastaFile = sys.argv[1]
renameFile = sys.argv[2]

curRealName = None
curPastaID = None
translate = {}
for i, line in enumerate(open(renameFile)):
    line = line.strip()
    if i % 3 == 0:
        curPastaID = line
    elif i % 3 == 1:
        curRealName = line
    else:
        translate[curPastaID] = curRealName

for header, seq in fastaRead(open(fastaFile)):
    # hacks for if we are using the badly-named original fasta.
    header = translate[header].replace("...", ".-.").replace(".", "_").replace("__", "_")
    fastaWrite(sys.stdout, header, seq)
