# Plots several comparable coalescence statistic XMLs.
require(plyr)
require(ggplot2)
require(XML)

# Data frame: genome1 genome2 identical early late identicalFraction earlyFraction lateFraction name
# Aggregate results are indicated by genome2="aggregate"
# Total aggregate results are genome1="aggregate", genome2="aggregate"
getCoalescenceResults <- function(xmlPath, name) {
    xmlFile <- xmlInternalTreeParse(xmlPath)
    # Genome-to-genome results
    genomeResults <- data.frame(t(sapply(xpathApply(xmlFile, "//genomeCoalescenceTest/coalescenceResults", xmlAttrs), unlist)))

    # Aggregate results per-genome
    genomes <- unlist(xpathApply(xmlFile, "//genomeCoalescenceTest", xmlGetAttr, "genome"))
    genomeAggregateResults <- data.frame(t(sapply(xpathApply(xmlFile, "//genomeCoalescenceTest/aggregateCoalescenceResults", xmlAttrs), unlist)))
    genomeAggregateResults$genome1 <- genomes
    genomeAggregateResults$genome2 <- rep("aggregate", length(genomeAggregateResults$genome1))

    # Total aggregate results
    aggregateResults <- data.frame(t(sapply(xpathApply(xmlFile, "/coalescenceTest/aggregateCoalescenceResults", xmlAttrs), unlist)))
    aggregateResults$genome1 <- rep("aggregate", length(aggregateResults$identical))
    aggregateResults$genome2 <- rep("aggregate", length(aggregateResults$identical))

    combined <- rbind(genomeResults, genomeAggregateResults, aggregateResults)
    combined$name <- rep(name, length(combined$genome1))
    combined <- transform(combined, identical=as.numeric(as.character(identical)), late=as.numeric(as.character(late)), early=as.numeric(as.character(early)),
              identicalFraction=as.numeric(as.character(identicalFraction)), lateFraction=as.numeric(as.character(lateFraction)), earlyFraction=as.numeric(as.character(earlyFraction)))
    return(combined)
}

df <- rbind(getCoalescenceResults("development.xml", "Development"), getCoalescenceResults("bestRecon_10trees_3out_nucLikelihood_1.0breakpoint.xml", "Breakpoint1.0"), getCoalescenceResults("bestRecon_10trees_3out_nucLikelihood_5.0breakpoint.xml", "Breakpoint5.0"), getCoalescenceResults("bestRecon_10trees_3out_nucLikelihood_noBreakpoint.xml", "Breakpoint0.0"))
pdf("test.pdf")
print(ggplot(subset(df, genome1 == "aggregate"), aes(y=identicalFraction, x=name)) + geom_bar(stat="identity") + theme_classic())
dev.off()
