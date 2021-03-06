#!/usr/bin/env python
"""Score the trees implied in a HAL file by sampling columns,
extracting surrounding sequence, realigning, and finally comparing
independently estimated trees vs. the induced trees from the HAL."""
from argparse import ArgumentParser
from jobTree.scriptTree.target import Target
from jobTree.scriptTree.stack import Stack
from sonLib.bioio import system, getTempFile, popenCatch
from sonLib.nxnewick import NXNewick
from collections import namedtuple, defaultdict, Counter
import math
import random
import sys

Coalescence = namedtuple('Coalescence', ['genome1', 'seq1', 'pos1', 'genome2', 'seq2', 'pos2', 'mrca'])

def getMRCA(tree, id1, id2):
    """Return the MRCA of two nodes in an NXTree."""
    candidates1 = set([id1])
    curNode = id1
    while tree.hasParent(curNode):
        curNode = tree.getParent(curNode)
        candidates1.add(curNode)
    curNode = id2
    if curNode in candidates1:
        return curNode
    while tree.hasParent(curNode):
        curNode = tree.getParent(curNode)
        if curNode in candidates1:
            return curNode
    raise RuntimeError("No MRCA found for nodes %d and %d" % (id1, id2))

def getNameToIdDict(tree):
    """Get a mapping from name to id for an nxtree. Will be invalid if
    changes are made to the tree.
    """
    ret = {}
    for id in tree.postOrderTraversal():
        if tree.hasName(id):
            ret[tree.getName(id)] = id
    return ret

def getLeafNames(tree):
    ret = []
    for id in tree.postOrderTraversal():
        if tree.hasName(id) and len(tree.getChildren(id)) == 0:
            ret.append(tree.getName(id))
    return ret

ColumnEntry = namedtuple('ColumnEntry', ['genome', 'seq', 'pos'])

def parseColumnEntryFromString(s):
    """Gets a ColumnEntry struct from a leaf name in a tree."""
    return ColumnEntry(genome=s.split("|")[0].split(".")[0],
                       seq=".".join(s.split("|")[0].split(".")[1:]),
                       pos=int(s.split("|")[1]))

def sampleCoalescences(tree, maxCoalescences, sampleNonDuplicates, requiredPosition=None):
    def choose(n, k):
        return math.factorial(n) / (math.factorial(k) * math.factorial(n - k))
    nameToId = getNameToIdDict(tree)
    validNames = getLeafNames(tree)

    # Find all duplicated genomes in a column
    numGenomeAppearances = Counter(parseColumnEntryFromString(i).genome for i in validNames)
    duplicatedGenomes = set(k for k, v in numGenomeAppearances.items() if v > 1)

    # Sample up to maxCoalescences pairs. This could be very slow if
    # maxCoalescences is high.
    pairs = set()
    visited = set()
    if requiredPosition is None:
        numPairs = choose(len(validNames), 2)
    else:
        numPairs = len(validNames) - 1
    numSamples = 0
    while numSamples < maxCoalescences and len(visited) != numPairs:
        pair = random.sample(validNames, 2)
        if requiredPosition is not None and pair[1] != "%s.%s|%s" % requiredPosition:
            pair[0] = "%s.%s|%s" % requiredPosition
        pair.sort()
        pair = (pair[0], pair[1])
        visited.add(pair)
        duplicated = parseColumnEntryFromString(pair[0]).genome in duplicatedGenomes \
                     and parseColumnEntryFromString(pair[1]).genome in duplicatedGenomes
        if pair not in pairs and (sampleNonDuplicates or duplicated):
            pairs.add(pair)
            numSamples += 1

    # Create coalescences out of the pairs
    coalescences = []
    for pair in pairs:
        mrca = tree.getName(getMRCA(tree, nameToId[pair[0]], nameToId[pair[1]]))
        # Relies on the sequences being named by
        # getRegionAroundSampledColumn, i.e. genome.seq|pos
        entry1 = parseColumnEntryFromString(pair[0])
        entry2 = parseColumnEntryFromString(pair[1])
        coalescence = Coalescence(genome1=entry1.genome, seq1=entry1.seq,
                                  pos1=entry1.pos,
                                  genome2=entry2.genome, seq2=entry2.seq,
                                  pos2=entry2.pos,
                                  mrca=mrca)
        coalescences.append(coalescence)
    return coalescences

def matchCoalescences(tree, inputCoalescences):
    """Find coalescences whose underlying pairs match the coalescences
    provided."""
    nameToId = getNameToIdDict(tree)
    coalescences = []
    for coalescence in inputCoalescences:
        # Relies on the sequences being named by
        # getRegionAroundSampledColumn, i.e. genome.seq|pos
        name1 = "%s.%s|%s" % (coalescence.genome1, coalescence.seq1, coalescence.pos1)
        name2 = "%s.%s|%s" % (coalescence.genome2, coalescence.seq2, coalescence.pos2)
        id1 = nameToId[name1]
        id2 = nameToId[name2]
        mrca = tree.getName(getMRCA(tree, id1, id2))
        coalescence = Coalescence(genome1=coalescence.genome1, seq1=coalescence.seq1, pos1=coalescence.pos1,
                                  genome2=coalescence.genome2, seq2=coalescence.seq2, pos2=coalescence.pos2,
                                  mrca=mrca)
        coalescences.append(coalescence)
    return coalescences

def getChromSizes(halPath, genome):
    """Get a dictionary of (chrom name):(chrom size) from a hal file."""
    output = popenCatch("halStats --chromSizes %s %s" % (genome, halPath))
    ret = {}
    for line in output.split("\n"):
        fields = line.split("\t")
        if len(fields) != 2:
            continue
        ret[fields[0]] = int(fields[1])
    return ret

def samplePosition(chromSizes):
    """Get a random (sequence, position) pair from a chromSizes dict. Very
    inefficient, but it shouldn't matter much.
    """
    genomeSize = sum(chromSizes.values()) # could pull this outside the function easily
    genomePos = random.randint(0, genomeSize - 1)
    curSize = 0
    for seq, size in chromSizes.items():
        if genomePos < curSize + size:
            return (seq, genomePos - curSize)
        curSize += size
    assert False

class Setup(Target):
    """Launch the sampling jobs and send the scores to the output
    phase."""
    def __init__(self, opts):
        Target.__init__(self)
        self.opts = opts

    def run(self):
        speciesTree = popenCatch("halStats --tree %s" % (self.opts.halFile)).strip()
        chromSizes = getChromSizes(self.opts.halFile, self.opts.refGenome)

        positions = []
        # For ensuring that a column isn't counted multiple times from
        # different reference positions.
        positionSet = set(positions)
        for i in xrange(self.opts.numSamples):
            # Have to sample the columns here since otherwise it can
            # be difficult to independently seed several RNGs
            pos = samplePosition(chromSizes)
            if pos not in positionSet:
                positions.append(pos)
                positionSet.add(pos)

        outputs = []
        for sliceStart in xrange(0, self.opts.numSamples,
                                 self.opts.samplesPerJob):
            slice = positions[sliceStart:sliceStart + self.opts.samplesPerJob]
            outputFile = getTempFile(rootDir=self.getGlobalTempDir())
            outputs.append(outputFile)
            self.addChildTarget(ScoreColumns(self.opts, slice,
                                             outputFile, speciesTree, positionSet))
        self.setFollowOnTarget(Summarize(self.opts, outputs, self.opts.outputFile, self.opts.writeMismatchesToFile))

class ScoreColumns(Target):
    """Get a column from the hal, realign the region surrounding the
    column, estimate a tree based on the realignment, then score the
    independently estimated tree against the one in the hal
    graph.
    """
    def __init__(self, opts, positions, outputFile, speciesTree, positionSet):
        Target.__init__(self)
        self.opts = opts
        self.positions = positions
        self.outputFile = outputFile
        self.speciesTree = speciesTree
        self.positionSet = positionSet

    def run(self):
        for position in self.positions:
            self.handleColumn(position)

    def handleColumn(self, position):
        # Get the column.
        fasta = popenCatch("getRegionAroundSampledColumn %s %s --refSequence %s --refPos %d" % (self.opts.halFile, self.opts.refGenome, position[0], position[1]))
        # Take out the tree (on the first line) in case the aligner is
        # picky (read: correct) about fasta parsing.
        fastaLines = fasta.split("\n")
        halTree = fastaLines[0][1:] # Skip '#' character.

        # Ensure that we only do each column once, by looking at the
        # ref genome's positions in the column. If this is not the
        # lowest position that was sampled, then we should stop to
        # avoid double-counting a column.
        headers = [l[1:] for l in fastaLines if len(l) > 0 and l[0] == '>']
        refGenomePoss = set((".".join(h.split("|")[0].split(".")[1:]), int(h.split("|")[-1])) for h in headers if h.split(".")[0] == self.opts.refGenome)
        if min(refGenomePoss.intersection(self.positionSet), key=lambda x: x[1]) != position:
            return

        # Check that the fasta actually has enough sequences to bother
        # with tree-building, and make sure it's duplicated if we want
        # only duplicated columns.
        seqNames = []
        for line in fastaLines:
            if len(line) != 0 and line[0] == '>':
                seqNames.append(line[1:])
        numGenomeAppearances = Counter(parseColumnEntryFromString(i).genome for i in seqNames)
        duplicatedGenomes = set(k for k, v in numGenomeAppearances.items() if v > 1)
        if len(seqNames) <= 3 or (not self.opts.nonDuplicated and len(duplicatedGenomes) == 0):
            return

        # Get rid of the initial comment line containing the newick tree.
        fasta = "\n".join(fastaLines[1:])

        # Align the region surrounding the column.
        alignInputPath = getTempFile(rootDir=self.getGlobalTempDir())
        open(alignInputPath, 'w').write(fasta)
        alignOutputPath = getTempFile(rootDir=self.getGlobalTempDir())
        alignCommand = self.opts.alignerCommand.replace("INPUT", alignInputPath).replace("OUTPUT", alignOutputPath)
        system(alignCommand)

        # Estimate a tree on the new alignment.
        treeOutputPath = getTempFile(rootDir=self.getGlobalTempDir())
        estimateCommand = self.opts.estimatorCommand.replace("INPUT", alignOutputPath).replace("OUTPUT", treeOutputPath)
        system(estimateCommand)
        estimatedTree = open(treeOutputPath).read().strip()

        # Reconcile against the species tree.

        # Make a spimap-esque "gene2species" file from the fact that
        # our sequences are labeled as genome.species|centerPos.
        gene2speciesPath = getTempFile(rootDir=self.getGlobalTempDir())
        gene2speciesHandle = open(gene2speciesPath, 'w')
        for line in fastaLines:
            if len(line) != 0 and line[0] == '>':
                header = line[1:]
                species = header.split(".")[0]
                gene2speciesHandle.write("%s\t%s\n" % (header, species))
        gene2speciesHandle.close()

        reconciled = popenCatch("reconcile '%s' '%s' '%s' 0 1" % (gene2speciesPath, estimatedTree, self.speciesTree))

        # Score the two trees
        self.reportCorrectCoalescences(position, halTree, reconciled)

    def reportCorrectCoalescences(self, position, halNewick, reconciledNewick):
        output = open(self.outputFile, 'a')
        hal = NXNewick().parseString(halNewick)
        reconciled = NXNewick().parseString(reconciledNewick)
        if self.opts.onlySelf:
            requiredPosition = ColumnEntry(self.opts.refGenome, position[0], position[1])
        else:
            requiredPosition = None
        halCoalescences = sampleCoalescences(hal, self.opts.coalescencesPerSample, self.opts.nonDuplicated, requiredPosition)
        reconciledCoalescences = matchCoalescences(reconciled, halCoalescences)
        assert(len(halCoalescences) == len(reconciledCoalescences))
        for halCoalescence, reconciledCoalescence in zip(halCoalescences, reconciledCoalescences):
            assert(halCoalescence.genome1 == reconciledCoalescence.genome1)
            assert(halCoalescence.seq1 == reconciledCoalescence.seq1)
            assert(halCoalescence.pos1 == reconciledCoalescence.pos1)
            assert(halCoalescence.genome2 == reconciledCoalescence.genome2)
            assert(halCoalescence.seq2 == reconciledCoalescence.seq2)
            assert(halCoalescence.pos2 == reconciledCoalescence.pos2)
            result = None
            speciesTree = NXNewick().parseString(self.speciesTree)
            nameToId = getNameToIdDict(speciesTree)
            # Have to get rid of the sequence/position information
            # in the hal MRCA
            halMrca = halCoalescence.mrca.split(".")[0]
            assert halMrca in nameToId
            reconciledId = nameToId[reconciledCoalescence.mrca]
            halId = nameToId[halMrca]
            id = getMRCA(speciesTree, halId, reconciledId)
            assert id == halId or id == reconciledId
            if reconciledId == halId:
                result = "identical"
            elif id == halId:
                # Late in hal relative to independent estimate
                result = "late"
            else:
                # Early in hal relative to independent estimate
                result = "early"
            if result != "identical" and self.opts.writeMismatchesToFile:
                output.write("mismatch\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (halCoalescence.genome1, halCoalescence.seq1, halCoalescence.pos1, halCoalescence.genome2, halCoalescence.seq2, halCoalescence.pos2, halMrca, reconciledCoalescence.mrca, result, halNewick, reconciledNewick))
            output.write("coalescence\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (halCoalescence.genome1, halCoalescence.seq1, halCoalescence.pos1, halCoalescence.genome2, halCoalescence.seq2, halCoalescence.pos2, result))

class CoalescenceResults:
    # Can't use namedtuple since tuples are immutable
    def __init__(self, identical=0, early=0, late=0):
        self.identical = 0
        self.early = 0
        self.late = 0

def defaultCoalescenceResults():
    return CoalescenceResults(0, 0, 0)

class Summarize(Target):
    """Merge the many output files into one large output file."""
    def __init__(self, opts, outputs, outputFile, mismatchPath):
        Target.__init__(self)
        self.opts = opts
        self.outputs = outputs
        self.outputFile = outputFile
        self.mismatchPath = mismatchPath

    def run(self):
        mismatchFile = None
        if self.mismatchPath is not None:
            mismatchFile = open(self.mismatchPath, 'w')
        results = defaultdict(lambda: defaultdict(defaultCoalescenceResults))
        results["aggregate"] = defaultCoalescenceResults()
        for output in self.outputs:
            for line in open(output):
                if line.strip() == "":
                    # blank line
                    continue
                fields = line.strip().split("\t")
                reportType = fields[0]
                fields = fields[1:]
                if reportType == "coalescence":
                    genome1 = fields[0]
                    seq1 = fields[1]
                    pos1 = fields[2]
                    genome2 = fields[3]
                    seq2 = fields[4]
                    pos2 = fields[5]
                    result = fields[6]
                    if result == "identical":
                        results["aggregate"].identical += 1
                        results[genome1]["aggregate"].identical += 1
                        results[genome2]["aggregate"].identical += 1
                        results[genome1][genome2].identical += 1
                        results[genome2][genome1].identical += 1
                    elif result == "early":
                        results["aggregate"].early += 1
                        results[genome1]["aggregate"].early += 1
                        results[genome2]["aggregate"].early += 1
                        results[genome1][genome2].early += 1
                        results[genome2][genome1].early += 1
                    else:
                        assert result == "late"
                        results["aggregate"].late += 1
                        results[genome1]["aggregate"].late += 1
                        results[genome2]["aggregate"].late += 1
                        results[genome2][genome1].late += 1
                        results[genome1][genome2].late += 1
                else:
                    print "\"%s\"" % reportType
                    print "%s" % line
                    assert reportType == "mismatch"
                    mismatchFile.write(line)
        with open(self.outputFile, 'w') as outputFile:
            outputFile.write('<coalescenceTest file="%s">\n' % (self.opts.halFile))
            self.printAggregateResults(outputFile, results["aggregate"])
            for genome1 in results.keys():
                if genome1 == "aggregate":
                    continue
                outputFile.write('<genomeCoalescenceTest genome="%s">\n' % (genome1))
                self.printAggregateResults(outputFile, results[genome1]["aggregate"])
                for genome2 in results[genome1]:
                    if genome2 == "aggregate":
                        continue
                    self.printGenomeResults(outputFile, genome1, genome2, results[genome1][genome2])
                outputFile.write('</genomeCoalescenceTest>\n')
            outputFile.write('</coalescenceTest>\n')

    def printAggregateResults(self, outputFile, results):
        total = results.identical + results.early + results.late
        if total == 0:
            # avoid a divide by zero. This just makes all the early,
            # late, identical fractions 0, which should be a sensible
            # response to getting no columns.
            total = 1
        outputFile.write('<aggregateCoalescenceResults identical="%d" early="%d" late="%d" identicalFraction="%f" earlyFraction="%f" lateFraction="%f" />\n' % (results.identical, results.early, results.late, float(results.identical)/total, float(results.early)/total, float(results.late)/total))

    def printGenomeResults(self, outputFile, genomeName1, genomeName2, results):
        total = results.identical + results.early + results.late
        if total == 0:
            # avoid a divide by zero. This just makes all the early,
            # late, identical fractions 0, which should be a sensible
            # response to getting no columns.
            total = 1
        outputFile.write('<coalescenceResults genome1="%s" genome2="%s" identical="%d" early="%d" late="%d" identicalFraction="%f" earlyFraction="%f" lateFraction="%f" />\n' % (genomeName1, genomeName2, results.identical, results.early, results.late, float(results.identical)/total, float(results.early)/total, float(results.late)/total))

if __name__ == '__main__':
    from scoreHalPhylogenies import * # required for jobTree
    parser = ArgumentParser(description=__doc__)
    Stack.addJobTreeOptions(parser)
    parser.add_argument('halFile', help='hal file')
    parser.add_argument('refGenome', help='reference genome')
    parser.add_argument('outputFile', help='output XML file')
    parser.add_argument('--numSamples', type=int,
                        help='Number of columns to sample',
                        default=50000)
    parser.add_argument('--nonDuplicated', action='store_true',
                        help='remove only-duplicated restriction on coalescence pairs',
                        default=False)
    parser.add_argument('--samplesPerJob', type=int,
                        help='Number of samples per jobTree job',
                        default=100)
    parser.add_argument('--coalescencesPerSample', type=int,
                        help='maximum number of coalescences to sample per column',
                        default=10)
    parser.add_argument('--onlySelf', default=False, action='store_true',
                        help='only sample coalescences including the sampled position')
    parser.add_argument('--width', type=int,
                        help='Width of region to extract around the'
                        ' sampled columns',  default=500)
    parser.add_argument('--alignerCommand', help='alignment command to run,'
                        ' where "INPUT" will be replaced with the input path'
                        ' and "OUTPUT" will be replaced with the output path',
                        default='mafft INPUT > OUTPUT')
    parser.add_argument('--estimatorCommand', help='tree-estimation command to run,'
                        ' where "INPUT" will be replaced with the input path'
                        ' and "OUTPUT" will be replaced with the output path',
                        default='fasttree -nt -gtr < INPUT > OUTPUT')
    parser.add_argument('--writeMismatchesToFile',
                        help="write trees to this file when at least one of "
                        "the sampled coalescences don't match")

    opts = parser.parse_args()
    Stack(Setup(opts)).startJobTree(opts)
