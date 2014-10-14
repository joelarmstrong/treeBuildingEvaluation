PROGRESSIVE_CACTUS_DIR = ~/progressiveCactus

.PHONY: all clean

all: bin/reconcile bin/getRegionAroundSampledColumn bin/guidedNeighborJoining

bin/reconcile: src/reconcile.c
	gcc -std=c99 -O3 -g -o $@ $< -I $(PROGRESSIVE_CACTUS_DIR)/submodules/sonLib/C/inc/ $(PROGRESSIVE_CACTUS_DIR)/submodules/sonLib/lib/*.a -lm -lstdc++

bin/getRegionAroundSampledColumn: src/getRegionAroundSampledColumn.cpp
	h5c++  -O3 -g -Wall -funroll-loops -DNDEBUG -I$(PROGRESSIVE_CACTUS_DIR)/submodules/sonLib/lib -I $(PROGRESSIVE_CACTUS_DIR)/submodules/hal/lib/ -o $@ $< $(PROGRESSIVE_CACTUS_DIR)/submodules/hal/lib/halLib.a $(PROGRESSIVE_CACTUS_DIR)/submodules/sonLib/lib/sonLib.a

bin/guidedNeighborJoining: src/guidedNeighborJoining.c
	gcc -std=c99 -O3 -g -o $@ $< -I $(PROGRESSIVE_CACTUS_DIR)/submodules/sonLib/C/inc/ $(PROGRESSIVE_CACTUS_DIR)/submodules/sonLib/lib/*.a -lm -lstdc++

clean:
	rm -f bin/*
