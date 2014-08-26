#!/usr/bin/env python
"""Build a MAF file out of a tree and a fasta file with sequences for
the leaves.

Usage: mafFromTreeAndFasta.py newickFile fastaFile > mafFile"""
import sys
from sonLib.bioio import fastaRead
from sonLib.nxnewick import NXNewick
from collections import defaultdict

seqLengths = {'chimpZnfCluster':2253077,
'gorillaZnfCluster':2337044,
'humanZnfCluster':2230928,
'orangZnfCluster':2367521,
'rhesusZnfCluster.fa':2174118}

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
    for header, seq in fastaRead(open(fastaPath)):
        fields = header.split('_')
        name = fields[0]
        start = fields[1]
        end = fields[2]
        strand = fields[3]
        sequences[header] = (seq, name, start, end, strand)
    
    # Print MAF, with sequence lines in post-order.
    print '##maf version=1 scoring=NA'
    print 'a tree="%s"' % (treeString)
    for nodeId in tree.postOrderTraversal():
        if not tree.isLeaf(nodeId):
            continue
        nodeName = tree.getName(nodeId)
        if nodeName not in sequences:
            raise RuntimeError("The tree has a node %s which was not found in the fasta file" % (nodeName))
        seq, name, start, end, strand = sequences[nodeName]
        alignedLen = lengthWithoutGaps(seq)
        seqLen = seqLengths[name]
        print 's %s %s %d %s %d %s' % (name, start, alignedLen, strand, seqLen, seq)
    # mafValidator wants an empty closing line(?)
    print ''
