PROGRESSIVE_CACTUS = ~/progressiveCactus/bin/runProgressiveCactus.sh
RUNS = $(wildcard run_*/)
HALS = $(addsuffix znfChr19.hal, $(RUNS))

all: $(HALS)

%/znfChr19.hal: %/cactus_progressive_config.xml %/znfChr19.txt
	rm -fr $(@D)/work
	$(PROGRESSIVE_CACTUS) --config $(@D)/cactus_progressive_config.xml $(@D)/znfChr19.txt $(@D)/work $@
