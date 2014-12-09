JOBTREE_ARGS = --maxThreads 45 --batchSystem parasol --parasolCommand '/cluster/home/jcarmstr/bin/parasol -host=ku'
PROGRESSIVE_CACTUS = ~/progressiveCactus/bin/runProgressiveCactus.sh
RUNS = $(wildcard run_*/)
HALS = $(addsuffix geneDesertish.hal, $(RUNS))
COALESCENCES = $(addsuffix humanCoalescences.xml, $(RUNS))
COVERAGES = $(addsuffix humanCoverage, $(RUNS))
HUBS = $(addsuffix hub, $(RUNS))
# because it's impossible to escape/quote spaces and commas in a makefile
SPACE :=
SPACE +=
COMMA := ,

all: $(HALS) $(COALESCENCES) $(COVERAGES) combined.pdf

hubs: $(HUBS)

%/geneDesertish.hal: %/cactus_progressive_config.xml %/geneDesertish.txt
	rm -fr $(@D)/work
	$(PROGRESSIVE_CACTUS) --config $(@D)/cactus_progressive_config.xml $(@D)/geneDesertish.txt $(@D)/work $@ --maxThreads 10 --stats

%/humanCoalescences.xml: %/geneDesertish.hal
	rm -fr $(@D)/jobTree
	PYTHONPATH=$(PYTHONPATH):../src PATH=$(PATH):../bin ../src/scoreHalPhylogenies.py --jobTree $(@D)/jobTree $(JOBTREE_ARGS) $< human $@
	rm -fr $(@D)/jobTree

%/humanCoverage: %/geneDesertish.hal
	halStats --inMemory --coverage human $< > $@

%/hub: %/geneDesertish.hal
	rm -fr $(@D)/hubJobTree
	hal2assemblyHub.py --lod --jobTree $(@D)/hubJobTree $< $@ --shortLabel geneDesertish$(@D) --longLabel geneDesertish$(@D)
	rm -fr $(@D)/hubJobTree

combined.pdf: $(COALESCENCES) $(COVERAGES)
	../src/plotCombinedStats.R $(subst $(SPACE),$(COMMA),$(COALESCENCES)) $(subst $(SPACE),$(COMMA),$(COVERAGES)) $(subst $(SPACE),$(COMMA),$(RUNS))
