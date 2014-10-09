JOBTREE_ARGS = --maxThreads 45 --batchSystem parasol --parasolCommand '/cluster/home/jcarmstr/bin/parasol -host=ku'
PROGRESSIVE_CACTUS = ~/progressiveCactus/bin/runProgressiveCactus.sh
RUNS = $(wildcard run_*/)
HALS = $(addsuffix znfChr19.hal, $(RUNS))
COALESCENCES = $(addsuffix humanCoalescences.xml, $(RUNS))
COVERAGES = $(addsuffix humanCoverage, $(RUNS))
HUBS = $(addsuffix hub, $(RUNS))
# because it's impossible to escape/quote spaces and commas in a makefile
SPACE :=
SPACE +=
COMMA := ,

all: $(HALS) $(COALESCENCES) $(COVERAGES) combined.pdf

hubs: $(HUBS)

%/znfChr19.hal: %/cactus_progressive_config.xml %/znfChr19.txt
	rm -fr $(@D)/work
	$(PROGRESSIVE_CACTUS) --config $(@D)/cactus_progressive_config.xml $(@D)/znfChr19.txt $(@D)/work $@ --maxThreads 10

%/humanCoalescences.xml: %/znfChr19.hal
	rm -fr $(@D)/jobTree
	PYTHONPATH=$(PYTHONPATH):.. PATH=$(PATH):.. ../scoreHalPhylogenies.py --jobTree $(@D)/jobTree $(JOBTREE_ARGS) $< human $@
	rm -fr $(@D)/jobTree

%/humanCoverage: %/znfChr19.hal
	halStats --inMemory --coverage human $< > $@

%/hub: %/znfChr19.hal
	rm -fr $(@D)/hubJobTree
	hal2assemblyHub.py --lod --bedDirs genes --jobTree $(@D)/hubJobTree $< $@ --shortLabel $(@D) --longLabel $(@D)
	rm -fr $(@D)/hubJobTree

combined.pdf: $(COALESCENCES) $(COVERAGES)
	../plotCombinedStats.R $(subst $(SPACE),$(COMMA),$(COALESCENCES)) $(subst $(SPACE),$(COMMA),$(COVERAGES)) $(subst $(SPACE),$(COMMA),$(RUNS))
