CACTUS_DIR = ~/cactus

.PHONY: all clean

all: bin/reconcile bin/getRegionAroundSampledColumn bin/guidedNeighborJoining

bin/reconcile: src/reconcile.c
	gcc -std=c99 -O3 -g -o $@ $< -I $(CACTUS_DIR)/submodules/sonLib/C/inc/ $(CACTUS_DIR)/submodules/sonLib/lib/*.a -lm -lstdc++

bin/getRegionAroundSampledColumn: src/getRegionAroundSampledColumn.cpp
	$(CACTUS_DIR)/submodules/hdf5/bin/h5c++  -O3 -g -Wall -funroll-loops -DNDEBUG -I $(CACTUS_DIR)/submodules/sonLib/C/inc/ -I $(CACTUS_DIR)/submodules/hal/api/inc -o $@ $< $(CACTUS_DIR)/submodules/hal/lib/libHal.a $(CACTUS_DIR)/submodules/sonLib/lib/sonLib.a

bin/guidedNeighborJoining: src/guidedNeighborJoining.c
	gcc -std=c99 -O3 -g -o $@ $< -I $(CACTUS_DIR)/submodules/sonLib/C/inc/ $(CACTUS_DIR)/submodules/sonLib/lib/*.a -lm -lstdc++

clean:
	rm -f bin/*
