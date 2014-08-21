#!/usr/bin/env python
# Filter chains by query coverage
import sys

if len(sys.argv) < 3:
    print "Usage: %s chain coveragePercent" % sys.argv[0]
chainFile = open(sys.argv[1])
minCoverage = float(sys.argv[2])

def outputChainIfPasses(chainLines, alignedLen, qSize, minCoverage):
    if 100*float(alignedLen)/qSize > minCoverage:
        print "\n".join(chainLines)
        print ''

curChain = []
qSize = 0
alignedLen = 0

for line in chainFile:
    line = line.strip()
    if line == '':
        continue
    if line[0] == "#":
        print line
        continue
    fields = line.split()
    if fields[0] == "chain":
        if len(curChain) > 0:
            outputChainIfPasses(curChain, alignedLen, qSize, minCoverage)
        # Chain header line
        curChain = []
        qSize = int(fields[8])
        alignedLen = 0
        curChain.append(line)
    else:
        # First field is length of ungapped block
        alignedLen += int(fields[0])
        curChain.append(line)

outputChainIfPasses(curChain, alignedLen, qSize, minCoverage)
