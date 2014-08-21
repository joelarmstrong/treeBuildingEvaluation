#include <stdlib.h>
#include <math.h>
#include "sonLib.h"
#include "stPinchPhylogeny.h"

static stHash *getLeafToSpecies(stTree *geneTree, stTree *speciesTree, stHash *leafNameToSpeciesName)
{
    stHash *ret = stHash_construct();
    stHashIterator *it = stHash_getIterator(leafNameToSpeciesName);
    char *leafName;
    while ((leafName = stHash_getNext(it)) != NULL) {
        stTree *gene = stTree_findChild(geneTree, leafName);
        assert(gene != NULL);
        char *speciesName = stHash_search(leafNameToSpeciesName, leafName);
        stTree *species = stTree_findChild(speciesTree, speciesName);
        stHash_insert(ret, gene, species);
    }
    return ret;
}

// Collapses internal nodes that have a single child.  Necessary since
// reconciliation is binary-only and we sometimes have a (somewhat
// sloppily-) induced tree.
static stTree *collapseUnnecessaryNodes(stTree *tree, double addToBranchLength)
{
    if (stTree_getChildNumber(tree) == 1) {
        addToBranchLength += stTree_getBranchLength(stTree_getChild(tree, 0));
        return collapseUnnecessaryNodes(stTree_getChild(tree, 0), addToBranchLength);
    } else {
        stTree *ret = stTree_cloneNode(tree);
        stTree_setBranchLength(ret, stTree_getBranchLength(tree) + addToBranchLength);
        for (int64_t i = 0; i < stTree_getChildNumber(tree); i++) {
            stTree_setParent(collapseUnnecessaryNodes(stTree_getChild(tree, i), 0), ret);
        }
        return ret;
    }
}

// Collapse internal nodes if they and their children have the same name.
static bool collapseIdenticalAncestors(stTree *tree)
{
    stTree *parent = stTree_getParent(tree);
    stTree_setBranchLength(tree, INFINITY);
    bool sameNameAsParent = false;

    if (parent != NULL) {
        sameNameAsParent = stString_eq(stTree_getLabel(tree), stTree_getLabel(parent));
    }

    for (int64_t i = 0; i < stTree_getChildNumber(tree); i++) {
        stTree *child = stTree_getChild(tree, i);
        if (collapseIdenticalAncestors(child)) {
            i = -1; // restart iteration.
        }
    }

    if (sameNameAsParent) {
        // remove this node.
        for (int64_t i = 0; i < stTree_getChildNumber(tree); i++) {
            stTree *child = stTree_getChild(tree, i);
            stTree_setParent(child, parent);
            i--;
        }
        stTree_setParent(tree, NULL);
        return true;
    }
    return false;
}

int main(int argc, char *argv[])
{
    int collapseIdenticalNodes = 0;
    if (argc < 4) {
        fprintf(stderr, "Usage: %s gene2species geneTree speciesTree [collapseIdenticalNodes=0]\n", argv[0]);
        return 1;
    }

    // read in the gene2species file.
    stHash *leafNameToSpeciesName = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, free, free);
    FILE *gene2speciesFile = fopen(argv[1], "r");
    char *line = malloc(sizeof(char));
    int64_t len = 1;
    while (benLine(&line, &len, gene2speciesFile) != -1) {
        stList *fields = stString_split(line);
        assert(stList_length(fields) == 2);
        stHash_insert(leafNameToSpeciesName, stString_copy(stList_get(fields, 0)), stString_copy(stList_get(fields, 1)));
        stList_destruct(fields);
    }

    stTree *geneTree = stTree_parseNewickString(argv[2]);
    geneTree = collapseUnnecessaryNodes(geneTree, 0.0);
    stTree *speciesTree = stTree_parseNewickString(argv[3]);
    stHash *leafToSpecies = getLeafToSpecies(geneTree, speciesTree, leafNameToSpeciesName);
    stPinchPhylogeny_reconcileAndLabelBinary(geneTree, speciesTree, leafToSpecies);

    if (argc > 4) {
        sscanf(argv[4], "%d", &collapseIdenticalNodes);
    }

    if (collapseIdenticalNodes) {
        collapseIdenticalAncestors(geneTree);
    }
    printf("%s\n", stTree_getNewickTreeString(geneTree));
}
