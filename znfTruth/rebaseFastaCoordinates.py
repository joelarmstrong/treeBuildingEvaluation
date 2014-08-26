#!/usr/bin/env python
import sys
from sonLib.bioio import fastaRead, fastaWrite
# Format: genome => sequence => subsequence, start of subsequence, end
# of subsequence
genomeMap = {"hg19":{"chr19":("humanZnfCluster", 51927367, 54158296)}, "panTro4":{"chr19":("chimpZnfCluster", 56310088, 58563166)},
             "gorGor3":{"chr19":("gorillaZnfCluster", 48765939, 51102984)}, "ponAbe2":{"chr19":("orangZnfCluster", 53063439, 55430961)},
             "rheMac3":{"chr19":("rhesusZnfCluster", 57314791, 59488909)}}

fasta = sys.argv[1]
renameFile = open(sys.argv[2], 'w')

for header, seq in fastaRead(open(fasta)):
    oldHeader = header
    fields = header.split("_")
    if len(fields) != 5:
        # some sequences have _'s in them (chrX_random_N)
        fields[0] = "_".join(fields[:len(fields)-4])
        fields = [field for i, field in enumerate(fields) if i == 0 or i > len(fields) - 5]
        assert len(fields) == 5
    chr = fields[0]
    start = int(fields[1])
    end = int(fields[2])
    strand = fields[3]
    genome = fields[4]
    if genome in genomeMap:
        if chr in genomeMap[genome]:
            subseq = genomeMap[genome][chr][0]
            subseqStart = genomeMap[genome][chr][1]
            subseqEnd = genomeMap[genome][chr][2]
            if start < subseqStart or end >= subseqEnd:
                continue
            assert start >= subseqStart
            assert end < subseqEnd
            newStart = None
            newEnd = None
            if strand == '+':
                newStart = start - subseqStart
                newEnd = end - subseqStart
            else:
                assert strand == '-'
                subseqSize = subseqStart - subseqEnd
                newStart = subseqEnd - end
                newEnd = subseqEnd - start
            header = "%s_%s_%s_%s" % (subseq, newStart, newEnd, strand)
            fastaWrite(sys.stdout, header, seq)
            renameFile.write("%s\n%s\n\n" % (oldHeader, header))
