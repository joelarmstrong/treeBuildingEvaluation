PROGRESSIVE_CACTUS_DIR = ~/progressiveCactus

all: reconcile getRegionAroundSampledColumn

reconcile: reconcile.c
	gcc -std=c99 -O0 -g -o reconcile reconcile.c -I $(PROGRESSIVE_CACTUS_DIR)/submodules/sonLib/C/inc/ -I  $(PROGRESSIVE_CACTUS_DIR)/submodules/pinchesAndCacti/inc $(PROGRESSIVE_CACTUS_DIR)/submodules/sonLib/lib/stPinchesAndCacti.a $(PROGRESSIVE_CACTUS_DIR)/submodules/sonLib/lib/*.a -lm -lstdc++

getRegionAroundSampledColumn: getRegionAroundSampledColumn.cpp
	h5c++  -O3 -g -Wall -funroll-loops -DNDEBUG -I$(PROGRESSIVE_CACTUS_DIR)/submodules/sonLib/lib -I $(PROGRESSIVE_CACTUS_DIR)/submodules/hal/lib/ -o getRegionAroundSampledColumn getRegionAroundSampledColumn.cpp $(PROGRESSIVE_CACTUS_DIR)/submodules/hal/lib/halLib.a $(PROGRESSIVE_CACTUS_DIR)/submodules/sonLib/lib/sonLib.a
