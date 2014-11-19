JOBTREE_ARGS = --maxThreads 45 --batchSystem parasol --parasolCommand '/cluster/home/jcarmstr/bin/parasol -host=ku'
PROGRESSIVE_CACTUS = ~/progressiveCactus/bin/runProgressiveCactus.sh
RUNS = $(wildcard run_*/)
HALS = $(addsuffix 6q14.2.hal, $(RUNS))
COALESCENCES = $(addsuffix humanCoalescences.xml, $(RUNS))
COVERAGES = $(addsuffix humanCoverage, $(RUNS))
HUBS = $(addsuffix hub, $(RUNS))
# because it's impossible to escape/quote spaces and commas in a makefile
SPACE :=
SPACE +=
COMMA := ,

all: $(HALS) $(COALESCENCES) $(COVERAGES) combined.pdf

hubs: $(HUBS)

%/6q14.2.hal: %/cactus_progressive_config.xml %/6q14.2.txt
	rm -fr $(@D)/work
	$(PROGRESSIVE_CACTUS) --config $(@D)/cactus_progressive_config.xml $(@D)/6q14.2.txt $(@D)/work $@ --maxThreads 10

%/humanCoalescences.xml: %/6q14.2.hal
	rm -fr $(@D)/jobTree
	PYTHONPATH=$(PYTHONPATH):../src PATH=$(PATH):../bin ../src/scoreHalPhylogenies.py --jobTree $(@D)/jobTree $(JOBTREE_ARGS) $< human $@
	rm -fr $(@D)/jobTree

%/humanCoverage: %/6q14.2.hal
	halStats --inMemory --coverage human $< > $@

%/hub: %/6q14.2.hal
	rm -fr $(@D)/hubJobTree
	hal2assemblyHub.py --lod --bedDirs genes --jobTree $(@D)/hubJobTree $< $@ --shortLabel 6q14.2$(@D) --longLabel 6q14.2$(@D)
	rm -fr $(@D)/hubJobTree

combined.pdf: $(COALESCENCES) $(COVERAGES)
	../src/plotCombinedStats.R $(subst $(SPACE),$(COMMA),$(COALESCENCES)) $(subst $(SPACE),$(COMMA),$(COVERAGES)) $(subst $(SPACE),$(COMMA),$(RUNS))
