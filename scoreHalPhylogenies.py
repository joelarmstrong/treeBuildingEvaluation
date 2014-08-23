#!/usr/bin/env python
"""Score the trees implied in a HAL file by sampling columns,
extracting surrounding sequence, realigning, and finally comparing
independently estimated trees vs. the induced trees from the HAL."""
from argparse import ArgumentParser
from jobTree.scriptTree.target import Target
from jobTree.scriptTree.stack import Stack
from sonLib.bioio import system, getTempFile, popenCatch

class Setup(Target):
    """Launch the sampling jobs and send the scores to the output
    phase."""
    def __init__(self, opts):
        Target.__init__(self)
        self.opts = opts

    def run(self):
        speciesTree = popenCatch("halStats --tree %s" % (self.opts.halFile)).strip()

        outputFile = getTempFile(rootDir=self.getGlobalTempDir())
        for i in xrange(self.opts.numSamples):
            self.addChildTarget(SampleAndScoreColumn(self.opts, outputFile,
                                                     speciesTree))
#        self.setFollowOnTarget(Output(self.opts, outputFile))

class SampleAndScoreColumn(Target):
    """Sample a column from the hal, realign the region surrounding the
    column, estimate a tree based on the realignment, then score the
    independently estimated tree against the one implied by the hal
    graph.
    """
    def __init__(self, opts, outputFile, speciesTree):
        Target.__init__(self)
        self.opts = opts
        self.outputFile = outputFile
        self.speciesTree = speciesTree

    def run(self):
        invalidColumn = True
        fasta = None
        while invalidColumn:
            invalidColumn = False
            # Sample a column.
            fasta = popenCatch("./getRegionAroundSampledColumn %s %s" % (self.opts.halFile, self.opts.refGenome))
            # Take out the tree (on the first line) in case the aligner is
            # picky (read: correct) about fasta parsing.
            fastaLines = fasta.split("\n")
            halTree = fastaLines[0][1:] # Skip '#' character.
            numSeqs = 0
            for line in fastaLines:
                if len(line) != 0 and line[0] == '>':
                    numSeqs += 1
            if numSeqs < 3:
                invalidColumn = True
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
        for line in fasta.split("\n"):
            if len(line) != 0 and line[0] == '>':
                header = line[1:]
                species = header.split(".")[0]
                gene2speciesHandle.write("%s\t%s\n" % (header, species))
        gene2speciesHandle.close()

        reconciled = popenCatch("./reconcile '%s' '%s' '%s' 1 1" % (gene2speciesPath, estimatedTree, self.speciesTree))
        self.logToMaster(reconciled)

class Output(Target):
    pass

if __name__ == '__main__':
    from scoreHalPhylogenies import * # required for jobTree
    parser = ArgumentParser(description=__doc__)
    Stack.addJobTreeOptions(parser)
    parser.add_argument('halFile', help='hal file')
    parser.add_argument('refGenome', help='reference genome')
    parser.add_argument('--numSamples', type=int,
                        help='Number of columns to sample',
                        default=100)
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

    opts = parser.parse_args()
    Stack(Setup(opts)).startJobTree(opts)
