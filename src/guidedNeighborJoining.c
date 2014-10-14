// Run guided neighbor-joining on a input (aligned!) FASTA file using
// the given join costs and output a tree.
#include <getopt.h>
#include "sonLib.h"
#include "bioioC.h"

// globals (necessary for fastaReadToFunction)
stHash *globalGenomeToMatrixIndex;
stHash *globalMatrixIndexToHeader;
int64_t globalSeqNum;

// Used with fastaReadToFunction to fill globalGenomeToMatrixIndex
// with stLists of stIntTuples representing genes' eventual position
// in the matrix index, and globalMatrixIndexToHeader keyed by
// stIntTuples and populated by their corresponding headers.
static void populateHashes(const char *header, const char *seq, int64_t length)
{
    // The genome is the first "."-separated token.
    stList *tokens = stString_splitByString(header, ".");
    if (stList_length(tokens) < 2) {
        st_errAbort("Header %s is not in UCSC genome.chr format.", header);
    }
    char *genome = stList_get(tokens, 0);

    if (stHash_search(globalGenomeToMatrixIndex, genome) == NULL) {
        // Initialize genome-keyed value with an empty list.
        stHash_insert(globalGenomeToMatrixIndex, stString_copy(genome), stList_construct3(0, (void (*)(void *)) stIntTuple_destruct));
    }
    stList *list = stHash_search(globalGenomeToMatrixIndex, genome);
    assert(list != NULL);
    stList_append(list, stIntTuple_construct1(globalSeqNum));
    stList_destruct(tokens);

    stHash_insert(globalMatrixIndexToHeader, stIntTuple_construct1(globalSeqNum), stString_copy(header));
    globalSeqNum++;
}

// Relabel a tree labeled by matrix indices in-place.
static void relabelTree(stTree *tree)
{
    for (int64_t i = 0; i < stTree_getChildNumber(tree); i++) {
        relabelTree(stTree_getChild(tree, i));
    }
    if (stTree_getLabel(tree) != 0) {
        assert(stTree_getChildNumber(tree) == 0);
        int64_t matrixIndex;
        assert(sscanf(stTree_getLabel(tree), "%" PRIi64, &matrixIndex) > 0);
        stIntTuple *queryIndex = stIntTuple_construct1(matrixIndex);
        char *header = stHash_search(globalMatrixIndexToHeader, queryIndex);
        assert(header != NULL);
        stTree_setLabel(tree, header);
        stIntTuple_destruct(queryIndex);
    }
}

// Check if either i or j are gap or wildcard characters.
static bool areNucsWildcards(char i, char j)
{
    i = tolower(i);
    j = tolower(j);
    if (i == 'n' || i == '-') {
        return true;
    }
    if (j == 'n' || j == '-') {
        return true;
    }
    return false;
}

// Check if i and j are equal despite masking.
static bool areNucsEqual(char i, char j)
{
    i = tolower(i);
    j = tolower(j);
    assert(i == 'a' || i == 'c' || i == 't' || i == 'g');
    assert(j == 'a' || j == 'c' || j == 't' || j == 'g');
    if (i == j) {
        return true;
    }
    return false;
}

// Add the number of similarities and differences between each of the
// non-gap column entries to the matrix.
static void addColumnInfoToMatrix(char *column, int64_t numEntries,
                                  stMatrix *matrix)
{
    for (int64_t i = 0; i < numEntries; i++) {
        for (int64_t j = i + 1; j < numEntries; j++) {
            // Check if either position i or j are Ns or -s.
            if (!areNucsWildcards(column[i], column[j])) {
                if (areNucsEqual(column[i], column[j])) {
                    // Mark this as a similarity.
                    *stMatrix_getCell(matrix, i, j) += 1.0;
                } else {
                    // Mark this as a difference.
                    *stMatrix_getCell(matrix, j, i) += 1.0;
                }
            }
        }
    }
}

static void usage(const char *progName)
{
    fprintf(stderr, "Usage: %s fastaFile speciesNewickTree\n", progName);
    fprintf(stderr, "Run neighbor-joining guided by a species tree on an "
            "aligned FASTA file of more than 3 DNA sequences. The sequences "
            "must be named in UCSC genome.chr format, and the species tree "
            "should be labeled by the genome names.\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, " --dupCost <float> : join cost per dup (default 0)\n");
    fprintf(stderr, " --lossCost <float> : join cost per loss (default 0)\n");
}

int main(int argc, char *argv[])
{
    // Parse options.
    double dupCost = 0.0;
    double lossCost = 0.0;
    struct option opts[] = {{"dupCost", required_argument, NULL, 'd'},
                            {"lossCost", required_argument, NULL, 'l'},
                            {NULL, 0, NULL, 0}};
    int optval = -1;
    while ((optval = getopt_long(argc, argv, "", opts, NULL)) != -1) {
        int valsParsed = 0;
        switch (optval) {
        case 'd':
            valsParsed = sscanf(optarg, "%lf", &dupCost);
            assert(valsParsed == 1);
            break;
        case 'l':
            valsParsed = sscanf(optarg, "%lf", &lossCost);
            assert(valsParsed == 1);
            break;
        case '?':
        default:
            usage(argv[0]);
            return 1;
        }
    }

    if (argc - optind < 2) {
        usage(argv[0]);
        return 1;
    }

    const char *fastaPath = argv[optind];
    const char *speciesNewick = argv[optind + 1];

    // Get columns from the pre-aligned fasta file.
    globalGenomeToMatrixIndex = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, free, (void (*)(void *)) stList_destruct);
    globalMatrixIndexToHeader = stHash_construct3((uint64_t (*)(const void *)) stIntTuple_hashKey, (int (*)(const void *, const void *)) stIntTuple_equalsFn, (void (*)(void *)) stIntTuple_destruct, free);
    globalSeqNum = 0;
    FILE *fastaFile = fopen(fastaPath, "r");
    fastaReadToFunction(fastaFile, populateHashes);
    fclose(fastaFile);
    struct CharColumnAlignment *columns = multiFastaRead((char *) fastaPath);

    // Create similarity matrix from columns.
    stMatrix *matrix = stMatrix_construct(columns->seqNo, columns->seqNo);
    for (int64_t i = 0; i < columns->columnNo; i++) {
        char *column = charColumnAlignment_getColumn(columns, i);
        addColumnInfoToMatrix(column, columns->seqNo, matrix);
    }

    // Parse species tree and create join cost matrix.
    stTree *speciesTree = stTree_parseNewickString(speciesNewick);
    stHash *speciesToJoinCostIndex = stHash_construct2(NULL, (void (*)(void *)) stIntTuple_destruct);
    stMatrix *joinCosts = stPhylogeny_computeJoinCosts(speciesTree, speciesToJoinCostIndex, dupCost, lossCost);

    // Map matrix indices to join cost indices.
    stHash *matrixIndexToJoinCostIndex = stHash_construct3((uint64_t (*)(const void *)) stIntTuple_hashKey, (int (*)(const void *, const void *)) stIntTuple_equalsFn, (void (*)(void *)) stIntTuple_destruct, (void (*)(void *)) stIntTuple_destruct);
    stHashIterator *genomeIt = stHash_getIterator(globalGenomeToMatrixIndex);
    char *curGenome;
    while ((curGenome = stHash_getNext(genomeIt)) != NULL) {
        stList *matrixIndices = stHash_search(globalGenomeToMatrixIndex, curGenome);
        assert(matrixIndices != NULL);
        stTree *species = stTree_findChild(speciesTree, curGenome);
        if (species == NULL) {
            st_errAbort("Node with name %s not found in species tree.", curGenome);
        }
        stIntTuple *joinCostIndex = stHash_search(speciesToJoinCostIndex, species);
        assert(joinCostIndex != NULL);
        for (int64_t i = 0; i < stList_length(matrixIndices); i++) {
            int64_t matrixIndex = stIntTuple_get(stList_get(matrixIndices, i), 0);
            stHash_insert(matrixIndexToJoinCostIndex, stIntTuple_construct1(matrixIndex), stIntTuple_construct1(stIntTuple_get(joinCostIndex, 0)));
        }
    }
    stHash_destructIterator(genomeIt);

    // Run the actual guided neighbor-joining.
    stTree *tree = stPhylogeny_guidedNeighborJoining(matrix, joinCosts, matrixIndexToJoinCostIndex, speciesToJoinCostIndex, speciesTree);
    relabelTree(tree);
    printf("%s\n", stTree_getNewickTreeString(tree));
    return 0;
}
