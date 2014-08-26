#!/usr/bin/env python
"""Build a MAF file out of a tree and a fasta file with sequences for
the leaves.

Usage: mafFromTreeAndFasta.py newickFile fastaFile > mafFile"""
import sys
from sonLib.bioio import fastaRead
from sonLib.nxnewick import NXNewick

def lengthWithoutGaps(seq):
    return len([i for i in seq if i != '-'])

if __name__ == '__main__':
    # Parse args
    if len(sys.argv) < 3:
        print __doc__
        sys.exit(1)

    newickPath = sys.argv[1]
    fastaPath = sys.argv[2]
    treeString = open(newickPath).read().split("\n")[0].strip()
    tree = NXNewick().parseString(treeString)
    
    sequences = {}
    for name, seq in fastaRead(open(fastaPath)):
        sequences[name] = seq
    
    # Print MAF, with sequence lines in post-order.
    print '##maf version=1 scoring=NA'
    print 'a tree="%s"' % (treeString)
    for nodeId in tree.postOrderTraversal():
        if not tree.isLeaf(nodeId):
            continue
        nodeName = tree.getName(nodeId)
        if nodeName not in sequences:
            raise RuntimeError("The tree has a node %s which was not found in the fasta file" % (nodeName))
        seq = sequences[nodeName]
        seqLen = lengthWithoutGaps(seq)
        print 's %s 0 %d + %d %s' % (nodeName, seqLen, seqLen, seq)
    # mafValidator wants an empty closing line(?)
    print ''
