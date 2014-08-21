#!/bin/bash
# Align and chain a ZNF gene to the target genome, and produce a BED
# file of likely ZNF genes.
set -o errexit
set -o nounset
set -o pipefail

# query file must be a 2bit
QUERYFILE=$1
QUERYSEQ=$2
TARGETGENOME=$3
MINCOVERAGE=$4
OUTPUTBED=$5
TEMPDIR=$(mktemp -d)
PATH=$PATH:/cluster/home/jcarmstr/bin/x86_64/:/cluster/home/jcarmstr/bin/

TARGETFILE=/hive/data/genomes/${TARGETGENOME}/${TARGETGENOME}.2bit

# Get chrom.sizes list from query
twoBitInfo ${QUERYFILE} ${TEMPDIR}/query.chrom.sizes

# Align
/cluster/home/jcarmstr/progressiveCactus/submodules/cactus/externalTools/lastz-distrib-1.03.54/src/lastz_32 --ambiguous=iupac ${TARGETFILE}[multiple] ${QUERYFILE}/${QUERYSEQ} --format=axt > ${TEMPDIR}/lastzLocal

# Chain
axtChain -linearGap=loose ${TEMPDIR}/lastzLocal ${TARGETFILE} ${QUERYFILE} ${TEMPDIR}/chain

# Filter by query coverage
/cluster/home/jcarmstr/znfChainTests/filterChainsByCoverage.py ${TEMPDIR}/chain ${MINCOVERAGE} > ${TEMPDIR}/filteredChains
# Filter so that the chains are relatively compact (have no large
# target gaps). Could also potentially filter using -tMaxSize.
chainFilter -tMaxGap=2000 ${TEMPDIR}/filteredChains > ${TEMPDIR}/compactChains

# Build bed file of all candidate genes in target
chainToPsl ${TEMPDIR}/compactChains /hive/data/genomes/${TARGETGENOME}/chrom.sizes ${TEMPDIR}/query.chrom.sizes ${TARGETFILE} ${QUERYFILE} ${TEMPDIR}/filtered.psl
pslToBed ${TEMPDIR}/filtered.psl ${TEMPDIR}/filtered.bed
bedtools bed12tobed6 -i ${TEMPDIR}/filtered.bed > ${TEMPDIR}/split.bed
bedtools sort -i ${TEMPDIR}/split.bed > ${TEMPDIR}/sorted.bed
bedtools merge -d 2000 -i ${TEMPDIR}/sorted.bed > ${OUTPUTBED}

# Clean up
rm -fr ${TEMPDIR}
