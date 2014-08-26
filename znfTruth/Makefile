# FIXME: make different targets so that this is actually useful
all:
	faSomeRecords -exclude pastaIteration1_origNames.fa pastaIteration1_sequenceBlacklist pastaIteration1_origNames.pruned.fa
	python rebaseFastaCoordinates.py pastaIteration1_origNames.pruned.fa origNamesToFinalNames > pastaIteration1_finalNames.fa
	python renameNewick.py origNamesToFinalNames pastaIteration1_origNames.nh > pastaIteration1_finalNames.nh
	python getSpimapGene2Species.py pastaIteration1_finalNames.fa >pastaIteration1_finalNames.smap
	gcc -std=c99 -O0 -g -o reconcile reconcile.c -I /cluster/home/jcarmstr/progressiveCactus/submodules/sonLib/C/inc/ -I  /cluster/home/jcarmstr/progressiveCactus/submodules/pinchesAndCacti/inc /cluster/home/jcarmstr/progressiveCactus/submodules/sonLib/lib/stPinchesAndCacti.a /cluster/home/jcarmstr/progressiveCactus/submodules/sonLib/lib/*.a -lm -lstdc++
	./reconcile pastaIteration1_finalNames.smap $(shell cat pastaIteration1_finalNames.nh) $(shell cat speciesTree.nh) 1 > pastaIteration1_finalNames.reconciled.nh
	 python mafFromTreeAndRebasedFasta.py pastaIteration1_finalNames.reconciled.nh pastaIteration1_finalNames.fa > pastaIteration1_toCompare.maf
