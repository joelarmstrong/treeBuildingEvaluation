// Samples a column and extracts a FASTA of the region surrounding
// every column entry.
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
  optionsParser->addOption("width", "width of the region around the sampled column to extract, default = 1000", 1000);
  return optionsParser;
}

int main(int argc, char *argv[])
{
  CLParserPtr optParser = initParser();
  string halPath, genomeName, refSequenceName;
  hal_index_t refPos = -1;
  hal_size_t width = 1000;
  try {
    optParser->parseOptions(argc, argv);
    halPath = optParser->getArgument<string>("halFile");
    genomeName = optParser->getArgument<string>("genome");
    refSequenceName = optParser->getOption<string>("refSequence");
    refPos = optParser->getOption<hal_index_t>("refPos");
    width = optParser->getOption<hal_size_t>("width");
  } catch (exception &e) {
    cerr << e.what() << endl;
    optParser->printUsage(cerr);
    return 1;
  }

  AlignmentConstPtr alignment = openHalAlignment(halPath, optParser);
  const Genome *genome = alignment->openGenome(genomeName);
  if (genome == NULL) {
    throw hal_exception("Genome " + genome->getName() + " not found in alignment");
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
  } else {
    refSequence = genome->getSequence(refSequenceName);
  }
  refPos += refSequence->getStartPosition();

  char *outputSeq = (char *) malloc((width*2 + 2) * sizeof(char));

  ColumnIteratorConstPtr colIt = genome->getColumnIterator(NULL, 0, refPos, NULL_INDEX, false, true);
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
      hal_size_t size = width * 2 + 1;
      if (size >= seq->getSequenceLength()) {
        size = seqEnd - seqStart;
      }

      // Be paranoid about the iterator still being reversed properly after we reposition it.
      bool reversed = dna->getReversed();
      dna->jumpTo(startPos);
      dna->setReversed(reversed);

      for(int64_t i = 0; i < (int64_t) size; i++, dna->toRight()) {
        outputSeq[i] = dna->getChar();
      }
      outputSeq[size] = '\0';

      stringstream header;
      header << genome->getName() << "." << seq->getName() << ":" << midpoint;
      fastaWrite(outputSeq, (char *) header.str().c_str(), stdout);
    }
  }
}
