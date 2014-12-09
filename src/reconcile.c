// Reconcile a gene tree with respect to a species tree, and,
// optionally, collapse duplication nodes to make a multifurcated
// tree.
#include <stdlib.h>
#include <math.h>
#include "bioioC.h"
#include "sonLib.h"

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
// sloppily-) induced tree with degree-2 nodes.
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

static bool isBinaryTree(stTree *tree)
{
    int64_t numChildren = stTree_getChildNumber(tree);
    if (numChildren != 0 && numChildren != 2) {
        return false;
    }
    for (int64_t i = 0; i < numChildren; i++) {
        stTree *child = stTree_getChild(tree, i);
        if (!isBinaryTree(child)) {
            return false;
        }
    }
    return true;
}

// Binarizes a tree that is at least binary (i.e. has at least 2
// children for every internal node).
static void arbitrarilyBinarize(stTree *tree)
{
    int64_t numChildren = stTree_getChildNumber(tree);
    if (isBinaryTree(tree)) {
        return;
    }
    assert(numChildren >= 2);

    for (int64_t i = 0; i < numChildren - 2; i++) {
        stTree *newNode = stTree_construct();
        stTree_setBranchLength(newNode, 0.0);
        stTree_setParent(newNode, tree);
        stTree_setParent(stTree_getChild(tree, 0), newNode);
        stTree_setParent(stTree_getChild(tree, 0), newNode);
    }
    assert(stTree_getChildNumber(tree) == 2);
    for (int64_t i = 0; i < 2; i++) {
        arbitrarilyBinarize(stTree_getChild(tree, i));
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
    int reRoot = 0;
    if (argc < 4) {
        fprintf(stderr, "Usage: %s gene2species geneTree speciesTree [collapseIdenticalNodes=0] [reRoot=0]\n", argv[0]);
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
    if (!isBinaryTree(geneTree)) {
        fprintf(stderr, "WARNING: arbitrarily binarizing multifurcated tree. This is "
                "fine if you are rerooting a trifurcated (unrooted) tree, but is bad "
                "otherwise.\n");
        arbitrarilyBinarize(geneTree);
        assert(isBinaryTree(geneTree));
    }
    stTree *speciesTree = stTree_parseNewickString(argv[3]);
    stHash *leafToSpecies = getLeafToSpecies(geneTree, speciesTree, leafNameToSpeciesName);

    if (argc > 4) {
        sscanf(argv[4], "%d", &collapseIdenticalNodes);
    }
    if (argc > 5) {
        sscanf(argv[5], "%d", &reRoot);
    }

    if (reRoot) {
        // If we actually cared about memory leaks this would be bad
        geneTree = stPhylogeny_rootByReconciliationAtMostBinary(geneTree, leafToSpecies);
        // Need to refresh the leafToSpecies map to correspond to the new tree
        leafToSpecies = getLeafToSpecies(geneTree, speciesTree, leafNameToSpeciesName);
    }

    // relabel the ancestors.
    stPhylogeny_reconcileAtMostBinary(geneTree, leafToSpecies, true);

    if (collapseIdenticalNodes) {
        collapseIdenticalAncestors(geneTree);
    }
    printf("%s\n", stTree_getNewickTreeString(geneTree));
}
