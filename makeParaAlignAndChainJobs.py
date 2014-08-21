#!/usr/bin/env python
from sonLib.bioio import system, popenCatch, getTempFile
import argparse
import os
# Quick script to output a bunch of jobs aligning and chaining from a
# 2bit to a target genome.

parser = argparse.ArgumentParser()
parser.add_argument("twoBit", help="2bit query file")
parser.add_argument("targetGenome", help="target genome name")
parser.add_argument("outputDir", help="Directory for output files")
parser.add_argument("--minCoverage", help="miniumum coverage %", default=80, type=float)

opts = parser.parse_args()

# Get list of seqs from the 2bit
seqInfo = popenCatch("twoBitInfo %s /dev/stdout" % opts.twoBit)
seqNames = map(lambda x: x[0], map(lambda y: y.split("\t"), seqInfo.split("\n")))
seqNames = filter(lambda x: x != "", seqNames)

system("mkdir -p %s" % opts.outputDir)

for seqName in seqNames:
    print "./alignAndChain.sh %s %s %s %f %s" % (opts.twoBit, seqName, opts.targetGenome, opts.minCoverage, os.path.join(opts.outputDir, seqName + ".bed"))
