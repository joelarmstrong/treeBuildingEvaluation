#!/usr/bin/env python
import sys
from sonLib.bioio import fastaRead
from sonLib.nxnewick import NXNewick
from sonLib.nxtree import NXTree
import networkx as nx

def induceTreeOnLeaves(nxtree, leaves):
    leaves = set(leaves)
    dg = nxtree.nxDg
    nodesToKeep = []
    for node in dg.nodes():
        succ = set([nxtree.getName(i) for i in nx.dfs_postorder_nodes(dg, node) if nxtree.hasName(i)])
        if len(succ.intersection(leaves)) != 0:
            nodesToKeep.append(node)
    return NXTree(dg.subgraph(nodesToKeep))

renameFile = open(sys.argv[1])
newickFile = open(sys.argv[2])

translate = {}
curPastaID = None
curRealName = None
for i, line in enumerate(renameFile):
    line = line.strip()
    if i % 3 == 0:
        curPastaID = line
    elif i % 3 == 1:
        curRealName = line
    else:
        translate[curPastaID] = curRealName.replace("...", ".-.").replace(".", "_").replace("__", "_")

s = newickFile.read()
for header1, header2 in translate.items():
    s = s.replace(header1, header2)

tree = NXNewick().parseString(s)

inducedTree = induceTreeOnLeaves(tree, translate.values())

print NXNewick().writeString(inducedTree)
