import sys
from sonLib.bioio import fastaRead
speciesMap = {'humanZnfCluster':'human', 'chimpZnfCluster':'chimp',
              'gorillaZnfCluster':'gorilla', 'rhesusZnfCluster.fa':'rhesus',
              'orangZnfCluster':'orang'}

for header, _ in fastaRead(open(sys.argv[1])):
    name = header.split("_")[0]
    print "%s\t%s" % (header, speciesMap[name])
