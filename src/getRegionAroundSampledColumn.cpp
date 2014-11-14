// Samples a column and extracts a FASTA of the region surrounding
// every column entry.
#include <time.h>
#include "hal.h"
#include "sonLib.h"
#include "bioioC.h"

using namespace std;
using namespace hal;

static CLParserPtr initParser()
{
    CLParserPtr optionsParser = hdf5CLParserInstance(true);
    optionsParser->addArgument("halFile", "target hal file");
    optionsParser->addArgument("genome", "reference genome");
    optionsParser->addOption("refSequence", "reference sequence name (by default, a random column will be sampled)", "");
    optionsParser->addOption("refPos", "position (only valid if also using --sequence)", -1);
    optionsParser->addOption("width", "width of the region around the sampled column to extract, default = 500", 500);
    optionsParser->addOptionFlag("lcaLabeling", "Label the ancestors with the MRCA of the child nodes instead of "
                                 "the ancestral genome and position (useful for comparing to a reconciled tree). "
                                 "NOTE: the ancestral sequences in the tree will not have their label changed.",
                                 false);
    return optionsParser;
}

// Adapted from the hal2maf block-tree building code.
// Could put this general version into the column iterator interface.

// Builds a "gene"-tree node and labels it properly.
static stTree *getTreeNode(SegmentIteratorConstPtr segIt)
{
    // Make sure the segment is sliced to only 1 base.
    assert(segIt->getStartPosition() == segIt->getEndPosition());
    stTree *ret = stTree_construct();
    const Genome *genome = segIt->getGenome();
    const Sequence *seq = genome->getSequenceBySite(segIt->getStartPosition());

    // make the sequence name safe for newick parsers (':'->'_')
    string seqName = seq->getName();
    hal_size_t i;
    while ((i = seqName.find(":")) != string::npos) {
        seqName[i] = '_';
    }

    stringstream ss;
    ss << segIt->getGenome()->getName() << "." << seqName << "|" << segIt->getStartPosition() - seq->getStartPosition();
    stTree_setLabel(ret, stString_copy(ss.str().c_str()));

    return ret;
}

// Recursive part of buildTree
// tree parameter represents node corresponding to the genome with
// bottom segment botIt
static void buildTreeR(BottomSegmentIteratorConstPtr botIt, stTree *tree)
{
    const Genome *genome = botIt->getGenome();

    // attach a node and recurse for each of this segment's children
    // (and paralogous segments)
    for (hal_size_t i = 0; i < botIt->getNumChildren(); i++) {
        if (botIt->hasChild(i)) {
            const Genome *child = genome->getChild(i);
            TopSegmentIteratorConstPtr topIt = child->getTopSegmentIterator();
            topIt->toChild(botIt, i);
            stTree *canonicalParalog = getTreeNode(topIt);
            stTree_setParent(canonicalParalog, tree);
            if (topIt->hasParseDown()) {
                BottomSegmentIteratorConstPtr childBotIt = child->getBottomSegmentIterator();
                childBotIt->toParseDown(topIt);
                buildTreeR(childBotIt, canonicalParalog);
            }
            // Traverse the paralogous segments cycle and add those segments as well
            if (topIt->hasNextParalogy()) {
                topIt->toNextParalogy();
                while(!topIt->isCanonicalParalog()) {
                    stTree *paralog = getTreeNode(topIt);
                    stTree_setParent(paralog, tree);
                    if(topIt->hasParseDown()) {
                        BottomSegmentIteratorConstPtr childBotIt = child->getBottomSegmentIterator();
                        childBotIt->toParseDown(topIt);
                        buildTreeR(childBotIt, paralog);
                    }
                    topIt->toNextParalogy();
                }
            }
        }
    }

    if (genome->getNumChildren() != 0 && stTree_getChildNumber(tree) == 0) {
        // Ancestral insertion. Ignore it.
        stTree_setParent(tree, NULL);
        stTree_destruct(tree);
    }
}

// Build a gene-tree from a column iterator.
static stTree *buildTree(ColumnIteratorConstPtr colIt)
{
    // Get any base from the column to begin building the tree
    const ColumnIterator::ColumnMap *colMap = colIt->getColumnMap();
    ColumnIterator::ColumnMap::const_iterator colMapIt = colMap->begin();
    const Sequence *sequence = NULL;
    hal_index_t index = NULL_INDEX;
    while (colMapIt != colMap->end()) {
        if (!colMapIt->second->empty()) {
            // Found a non-empty column map entry, just take the index and
            // sequence of the first base found
            sequence = colMapIt->first;
            index = colMapIt->second->at(0)->getArrayIndex();
            break;
        }
        colMapIt++;
    }
    assert(sequence != NULL && index != NULL_INDEX);
    const Genome *genome = sequence->getGenome();

    // Get the bottom segment that is the common ancestor of all entries
    TopSegmentIteratorConstPtr topIt = genome->getTopSegmentIterator();
    BottomSegmentIteratorConstPtr botIt;
    if (genome->getNumTopSegments() == 0) {
        // The reference is the root genome.
        botIt = genome->getBottomSegmentIterator();
        botIt->toSite(index);
    } else {
        // Keep heading up the tree until we hit the root segment.
        topIt->toSite(index);
        while (topIt->hasParent()) {
            const Genome *parent = topIt->getGenome()->getParent();
            botIt = parent->getBottomSegmentIterator();
            botIt->toParent(topIt);
            if(parent->getParent() == NULL || !botIt->hasParseUp()) {
                // Reached root genome
                break;
            }
            topIt = parent->getTopSegmentIterator();
            topIt->toParseUp(botIt);
        }
    }

    stTree *tree = NULL;
    if(topIt->hasParent() == false && topIt->getGenome() == genome && genome->getNumBottomSegments() == 0) {
        // Handle insertions in leaves. botIt doesn't point anywhere since
        // there are no bottom segments.
        tree = getTreeNode(topIt);
    } else {
        tree = getTreeNode(botIt);
        buildTreeR(botIt, tree);
    }
    assert(tree != NULL);
    return tree;
}

static void relabelAncestorsToLCAOfChildren(stTree *tree, stTree *speciesTree) {
    if (stTree_getChildNumber(tree) == 0) {
        // Not relabeling leaves.
        return;
    }

    stTree *curMRCA = NULL;
    for (int64_t i = 0; i < stTree_getChildNumber(tree); i++) {
        relabelAncestorsToLCAOfChildren(stTree_getChild(tree, i), speciesTree);
        const char *childLabel = stTree_getLabel(stTree_getChild(tree, i));
        // Get the first .-separated token, i.e. the genome name.
        stList *tokens = stString_splitByString(childLabel, ".");
        char *childSpeciesLabel = (char *) stList_get(tokens, 0);
        stTree *childSpecies = stTree_findChild(speciesTree, childSpeciesLabel);
        assert(childSpecies != NULL);
        stList_destruct(tokens);
        if (curMRCA == NULL) {
            curMRCA = childSpecies;
        } else {
            curMRCA = stTree_getMRCA(curMRCA, childSpecies);
        }
    }
    stTree_setLabel(tree, stTree_getLabel(curMRCA));
}

int main(int argc, char *argv[])
{
    CLParserPtr optParser = initParser();
    string halPath, genomeName, refSequenceName;
    hal_index_t refPos = -1;
    hal_index_t width = 1000;
    bool lcaLabeling;
    try {
        optParser->parseOptions(argc, argv);
        halPath = optParser->getArgument<string>("halFile");
        genomeName = optParser->getArgument<string>("genome");
        refSequenceName = optParser->getOption<string>("refSequence");
        refPos = optParser->getOption<hal_index_t>("refPos");
        width = optParser->getOption<hal_index_t>("width");
        lcaLabeling = optParser->getFlag("lcaLabeling");
    } catch (exception &e) {
        cerr << e.what() << endl;
        optParser->printUsage(cerr);
        return 1;
    }

    st_randomSeed(time(NULL));

    AlignmentConstPtr alignment = openHalAlignmentReadOnly(halPath, optParser);
    const Genome *genome = alignment->openGenome(genomeName);
    if (genome == NULL) {
        throw hal_exception("Genome " + genomeName + " not found in alignment");
    }

    const Sequence *refSequence = NULL;
    if (refSequenceName.empty()) {
        // Sample a position from the entire genome.
        hal_size_t genomeLen = genome->getSequenceLength();
        refPos = st_randomInt64(0, genomeLen - 1);
        refSequence = genome->getSequenceBySite(refPos);
    } else if (refPos == -1) {
        refSequence = genome->getSequence(refSequenceName);
        hal_size_t sequenceLen = refSequence->getSequenceLength();
        refPos = st_randomInt64(0, sequenceLen - 1);
        refPos += refSequence->getStartPosition();
    } else {
        refSequence = genome->getSequence(refSequenceName);
        refPos += refSequence->getStartPosition();
    }

    char *outputSeq = (char *) malloc((width*2 + 2) * sizeof(char));

    ColumnIteratorConstPtr colIt = genome->getColumnIterator(NULL, 0, refPos, NULL_INDEX, false, true);

    stTree *colTree = buildTree(colIt);
    if (lcaLabeling) {
        stTree *speciesTree = stTree_parseNewickString(alignment->getNewickTree().c_str());
        relabelAncestorsToLCAOfChildren(colTree, speciesTree);
    }

    // Print out the tree underlying this column as a FASTA comment.
    cout << "#" << stTree_getNewickTreeString(colTree) << endl;

    const ColumnIterator::ColumnMap *cols = colIt->getColumnMap();
    ColumnIterator::ColumnMap::const_iterator colMapIt;
    for (colMapIt = cols->begin(); colMapIt != cols->end(); colMapIt++) {
        if (colMapIt->second->empty()) {
            // The column map can contain empty entries.
            continue;
        }

        const Sequence *seq = colMapIt->first;
        hal_index_t seqStart = seq->getStartPosition();
        hal_index_t seqEnd = seqStart + seq->getSequenceLength(); // exclusive.
        for (ColumnIterator::DNASet::const_iterator dnaIt = colMapIt->second->begin(); dnaIt != colMapIt->second->end(); dnaIt++) {
            DNAIteratorConstPtr dna = *dnaIt;
            hal_index_t midpoint = dna->getArrayIndex();
            hal_index_t startPos = (dna->getReversed()) ? midpoint + width : midpoint - width;
            if (startPos < seqStart) {
                startPos = seqStart;
            } else if (startPos >= seqEnd) {
                startPos = seqEnd - 1;
            }
            hal_index_t size = width * 2 + 1;
            if ((hal_size_t) size >= seq->getSequenceLength()) {
                size = seqEnd - seqStart;
            }

            if (dna->getReversed() == false && startPos + size >= seqEnd) {
                size = seqEnd - startPos - 1;
            } else if(dna->getReversed() == true && startPos - size < seqStart) {
                size = startPos - seqStart - 1;
            }

            // Be paranoid about the iterator still being reversed properly after we reposition it.
            bool reversed = dna->getReversed();
            dna->jumpTo(startPos);
            dna->setReversed(reversed);

            for(int64_t i = 0; i < (int64_t) size; i++, dna->toRight()) {
                outputSeq[i] = dna->getChar();
            }
            outputSeq[size] = '\0';

            // make the sequence name safe for newick parsers (':'->'_')
            string seqName = seq->getName();
            hal_size_t i;
            while ((i = seqName.find(":")) != string::npos) {
                seqName[i] = '_';
            }

            stringstream header;
            header << seq->getGenome()->getName() << "." << seqName << "|" << midpoint - seqStart;
            fastaWrite(outputSeq, (char *) header.str().c_str(), stdout);
        }
    }

    // Intentionally not dealing with memory leaks for this very
    // short-lived process.
}
