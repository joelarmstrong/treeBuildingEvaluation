PROGRESSIVE_CACTUS = ~/progressiveCactus/bin/runProgressiveCactus.sh
RUNS = $(wildcard run_*/)
HALS = $(addsuffix znfChr19.hal, $(RUNS))
COALESCENCES = $(addsuffix humanCoalescences.xml, $(RUNS))
COVERAGES = $(addsuffix humanCoverage, $(RUNS))

# because it's impossible to escape/quote spaces and commas in a makefile
SPACE :=
SPACE +=
COMMA := ,


all: $(HALS) $(COALESCENCES) $(COVERAGES) combined.pdf

%/znfChr19.hal: %/cactus_progressive_config.xml %/znfChr19.txt
	rm -fr $(@D)/work
	$(PROGRESSIVE_CACTUS) --config $(@D)/cactus_progressive_config.xml $(@D)/znfChr19.txt $(@D)/work $@ --maxThreads 10

%/humanCoalescences.xml: %/znfChr19.hal
	rm -fr $(@D)/jobTree
	PATH=$(PATH):.. ../scoreHalPhylogenies.py --jobTree $(@D)/jobTree --maxThreads 8 $< human $@

%/humanCoverage: %/znfChr19.hal
	halStats --inMemory --coverage human $< > $@

combined.pdf: $(COALESCENCES) $(COVERAGES)
	../plotCombinedStats.R $(subst $(SPACE),$(COMMA),$(COALESCENCES)) $(subst $(SPACE),$(COMMA),$(COVERAGES)) $(RUNS)
