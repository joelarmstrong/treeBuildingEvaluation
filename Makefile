PROGRESSIVE_CACTUS_DIR = ~/progressiveCactus

.PHONY: all clean

all: reconcile getRegionAroundSampledColumn guidedNeighborJoining

reconcile: reconcile.c
	gcc -std=c99 -O3 -g -o reconcile reconcile.c -I $(PROGRESSIVE_CACTUS_DIR)/submodules/sonLib/C/inc/ $(PROGRESSIVE_CACTUS_DIR)/submodules/sonLib/lib/*.a -lm -lstdc++

getRegionAroundSampledColumn: getRegionAroundSampledColumn.cpp
	h5c++  -O3 -g -Wall -funroll-loops -DNDEBUG -I$(PROGRESSIVE_CACTUS_DIR)/submodules/sonLib/lib -I $(PROGRESSIVE_CACTUS_DIR)/submodules/hal/lib/ -o getRegionAroundSampledColumn getRegionAroundSampledColumn.cpp $(PROGRESSIVE_CACTUS_DIR)/submodules/hal/lib/halLib.a $(PROGRESSIVE_CACTUS_DIR)/submodules/sonLib/lib/sonLib.a

guidedNeighborJoining: guidedNeighborJoining.c
	gcc -std=c99 -O3 -g -o guidedNeighborJoining guidedNeighborJoining.c -I $(PROGRESSIVE_CACTUS_DIR)/submodules/sonLib/C/inc/ $(PROGRESSIVE_CACTUS_DIR)/submodules/sonLib/lib/*.a -lm -lstdc++

clean:
	rm -f reconcile getRegionAroundSampledColumn guidedNeighborJoining
